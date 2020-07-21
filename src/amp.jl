"""

    create_bounding_mip(m::Optimizer; use_disc::Dict)

Set up a JuMP MILP bounding model base on variable domain partitioning information stored in `use_disc`.
By default, if `use_disc is` not provided, it will use `m.discretizations` store in the Alpine model.
The basic idea of this MILP bounding model is to use Tighten McCormick to convexify the original Non-convex region.
Among all presented partitions, the bounding model will choose one specific partition as the lower bound solution.
The more partitions there are, the better or finer bounding model relax the original MINLP while the more
efforts required to solve this MILP is required.

This function is implemented in the following manner:

    * [`amp_post_vars`](@ref): post original and lifted variables
    * [`amp_post_lifted_constraints`](@ref): post original and lifted constraints
    * [`amp_post_lifted_obj`](@ref): post original or lifted objective function
    * [`amp_post_tmc_mccormick`](@ref): post Tighten McCormick variables and constraints base on `discretization` information

More specifically, the Tightening McCormick used here can be generalized in the following mathematical formulation. Consider a nonlinear term
```math
\\begin{subequations}
\\begin{align}
   &\\widehat{x_{ij}} \\geq (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\widehat{x_{ij}} \\geq (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\widehat{x_{ij}} \\leq (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\widehat{x_{ij}} \\leq (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) \\\\
   & \\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i) \\geq x_{i} \\geq \\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i) \\\\
   & \\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) \\geq x_{j} \\geq \\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\sum \\hat{\\mathbf{y}_i} = 1, \\ \\ \\sum \\hat{\\mathbf{y}_j}_k = 1 \\\\
   &\\hat{\\mathbf{y}}_i \\in \\{0,1\\}, \\hat{\\mathbf{y}}_j \\in \\{0,1\\}
\\end{align}
\\end{subequations}
```

"""
function create_bounding_mip(m::Optimizer; use_disc=nothing)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc

    m.model_mip = Model(get_option(m, :mip_solver)) # Construct JuMP Model
    start_build = time()
    # ------- Model Construction ------ #
    amp_post_vars(m)                                                # Post original and lifted variables
    amp_post_lifted_constraints(m)                                  # Post lifted constraints
    amp_post_lifted_objective(m)                                    # Post objective
    amp_post_convexification(m, use_disc=discretization)            # Convexify problem
    # --------------------------------- #
    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, get_option(m, :timeout) - m.logs[:total_time])

    return
end

"""
    amp_post_convexification(m::Optimizer; kwargs...)

wrapper function to convexify the problem for a bounding model. This function talks to nonconvex_terms and convexification methods
to finish the last step required during the construction of bounding model.
"""
function amp_post_convexification(m::Optimizer; use_disc=nothing)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc

    for i in 1:length(get_option(m, :method_convexification))             # Additional user-defined convexification method
        eval(get_option(m, :method_convexification)[i])(m)
    end

    amp_post_mccormick(m, use_disc=discretization)          # handles all bi-linear and monomial convexificaitons
    amp_post_convhull(m, use_disc=discretization)           # convex hull representation

    is_fully_convexified(m) # Ensure if  all the non-linear terms are convexified

    return
end

function amp_post_vars(m::Optimizer; kwargs...)

    options = Dict(kwargs)

    if haskey(options, :use_disc)
        l_var = [options[:use_disc][i][1]   for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
        u_var = [options[:use_disc][i][end] for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
    else
        l_var = m.l_var_tight
        u_var = m.u_var_tight
    end

    @variable(m.model_mip, x[i=1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)])

    for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)
        # This is a tricky step, not enforcing category of lifted variables is able to improve performance
        (i <= m.num_var_orig) && setcategory(x[i], m.var_type_orig[i])
        # Changed to tight bound, if no bound tightening is performed, will be just .l_var_orig
        l_var[i] > -Inf && JuMP.set_lower_bound(x[i], l_var[i])
        # Changed to tight bound, if no bound tightening is performed, will be just .u_var_orig
        u_var[i] < Inf && JuMP.set_upper_bound(x[i], u_var[i])

        m.var_type[i] == :Int && error("Support for general integer problem is current limited...")
    end

    return
end


function amp_post_lifted_constraints(m::Optimizer)

    for i in 1:m.num_constr_orig
        if m.constr_structure[i] == :affine
            amp_post_affine_constraint(m.model_mip, m.bounding_constr_mip[i])
        elseif m.constr_structure[i] == :convex
            amp_post_convex_constraint(m.model_mip, m.bounding_constr_mip[i])
        else
            error("Unknown constr_structure type $(m.constr_structure[i])")
        end
    end

    for i in keys(m.linear_terms)
        amp_post_linear_lift_constraints(m.model_mip, m.linear_terms[i])
    end

    return
end

function amp_post_affine_constraint(model_mip::JuMP.Model, affine::Dict)

    if affine[:sense] == :(>=)
        @constraint(model_mip,
            sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) >= affine[:rhs])
    elseif affine[:sense] == :(<=)
        @constraint(model_mip,
            sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) <= affine[:rhs])
    elseif affine[:sense] == :(==)
        @constraint(model_mip,
            sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) == affine[:rhs])
    else
        error("Unkown sense.")
    end

    return
end

function amp_post_convex_constraint(model_mip::JuMP.Model, convex::Dict)

    !prod([i == 2 for i in convex[:powers]]) && error("No relaxation implementation for convex constraints $(convex[:expr])")

    if convex[:sense] == :(<=)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) <= convex[:rhs])
    elseif convex[:sense] == :(>=)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) >= convex[:rhs])
    else
        error("No equality constraints should be recognized as supported convex constriants")
    end

    return
end

function amp_post_linear_lift_constraints(model_mip::JuMP.Model, l::Dict)

    @assert l[:ref][:sign] == :+
    @constraint(model_mip, Variable(model_mip, l[:y_idx]) == sum(i[1]*Variable(model_mip, i[2]) for i in l[:ref][:coef_var]) + l[:ref][:scalar])
    return
end

function amp_post_lifted_objective(m::Optimizer)

#if isa(m.obj_expr_orig, Number)
if expr_isconst(m.obj_expr_orig)
    @objective(m.model_mip, m.sense_orig, eval(m.obj_expr_orig))
   elseif m.obj_structure == :affine
        @objective(m.model_mip, m.sense_orig, m.bounding_obj_mip[:rhs] + sum(m.bounding_obj_mip[:coefs][i]*Variable(m.model_mip, m.bounding_obj_mip[:vars][i].args[2]) for i in 1:m.bounding_obj_mip[:cnt]))
    elseif m.obj_structure == :convex
        # This works only when the original objective is convex quadratic.
        # Higher-order convex monomials need implementation of outer-approximation (check resolve_convex_constr in operators.jl)
        @objective(m.model_mip, m.sense_orig, m.bounding_obj_mip[:rhs] + sum(m.bounding_obj_mip[:coefs][i]*Variable(m.model_mip, m.bounding_obj_mip[:vars][i].args[2])^2 for i in 1:m.bounding_obj_mip[:cnt]))
    else
     @show m.obj_expr_orig
        error("Unknown structural obj type $(m.obj_structure)")
    end
    return
end

function add_partition(m::Optimizer; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization
    haskey(options, :use_solution) ? point_vec = options[:use_solution] : point_vec = m.best_bound_sol

    if isa(get_option(m, :disc_add_partition_method), Function)
        m.discretization = eval(get_option(m, :disc_add_partition_method))(m, use_disc=discretization, use_solution=point_vec)
    elseif get_option(m, :disc_add_partition_method) == "adaptive"
        m.discretization = add_adaptive_partition(m, use_disc=discretization, use_solution=point_vec)
    elseif get_option(m, :disc_add_partition_method) == "uniform"
        m.discretization = add_uniform_partition(m, use_disc=discretization)
    else
        error("Unknown input on how to add partitions.")
    end

    return
end

"""
    add_discretization(m::Optimizer; use_disc::Dict, use_solution::Vector)

Basic built-in method used to add a new partition on feasible domains of discretizing variables.
This method makes modification in discretization

Consider an original partition [0, 3, 7, 9], where LB/any solution is 4.
Use ^ as the new partition, "|" as the original partition

A case when discretize ratio = 4
| -------- | - ^ -- * -- ^ ---- | -------- |
0          3  3.5   4   4.5     7          9

A special case when discretize ratio = 2
| -------- | ---- * ---- ^ ---- | -------- |
0          3      4      5      7          9

There are two options for this function,

    * `use_disc(default=m.discretization)`:: to regulate which is the base to add new partitions on
    * `use_solution(default=m.best_bound_sol)`:: to regulate which solution to use when adding new partitions on

TODO: also need to document the special diverted cases when new partition touches both corners

This function can be accordingly modified by the user to change the behavior of the solver, and thus the convergence.
"""
function add_adaptive_partition(m::Optimizer;kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization
    haskey(options, :use_solution) ? point_vec = copy(options[:use_solution]) : point_vec = copy(m.best_bound_sol)
    haskey(options, :use_ratio) ? ratio = options[:use_ratio] : ratio = get_option(m, :disc_ratio)
    haskey(options, :branching) ? branching = options[:branching] : branching = false

    if length(point_vec) < m.num_var_orig + m.num_var_linear_mip + m.num_var_nonlinear_mip
        point_vec = resolve_lifted_var_value(m, point_vec)  # Update the solution vector for lifted variable
    end

    if branching
        discretization = deepcopy(discretization)
    end

    processed = Set{Int}()

    # ? Perform discretization based on type of nonlinear terms ? #
    for i in m.disc_vars
        point = point_vec[i]                # Original Variable
        λCnt = length(discretization[i])

        # @show i, point, discretization[i] # Debugger

        # Built-in method based-on variable type
        if m.var_type[i] == :Cont
            i in processed && continue
            point = correct_point(m, discretization[i], point, i)
            for j in 1:λCnt
                if point >= discretization[i][j] && point <= discretization[i][j+1]  # Locating the right location
                    radius = calculate_radius(discretization[i], j, ratio)
                    insert_partition(m, i, j, point, radius, discretization[i])
                    push!(processed, i)
                    break
                end
            end
        elseif m.var_type[i] == :Bin # This should never happen
            @warn "  Warning: Binary variable in m.disc_vars. Check out what is wrong..."
            continue  # No partition should be added to binary variable unless user specified
        elseif m.var_type[i] == :Int
            error("Support for general integer problem is current limited...")
        else
            error("Unexpected variable types during injecting partitions")
        end
    end

    return discretization
end

"""
    This function targets to address unexpected numerical issues when adding partitions in tight regions.
"""
function correct_point(m::Optimizer, partvec::Vector, point::Float64, var::Int)

    if point < partvec[1] - get_option(m, :tol) || point > partvec[end] + get_option(m, :tol)
        @warn "  Warning: VAR$(var) SOL=$(point) out of discretization [$(partvec[1]),$(partvec[end])]. Taking middle point..."
        return 0.5*(partvec[1] + partvec[end]) # Should choose the longest range
    end

    isapprox(point, partvec[1];atol=get_option(m, :tol)) && return partvec[1]
    isapprox(point, partvec[end];atol=get_option(m, :tol)) && return partvec[end]

    return point
end

function calculate_radius(partvec::Vector, part::Int, ratio::Any)

    lb_local = partvec[part]
    ub_local = partvec[part+1]

    distance = ub_local - lb_local
    if isa(ratio, Float64) || isa(ratio, Int)
        radius = distance / ratio
    elseif isa(ratio, Function)
        radius = distance / ratio(m)
    else
        error("Undetermined discretization ratio")
    end

    return radius
end

function insert_partition(m::Optimizer, var::Int, partidx::Int, point::Number, radius::Float64, partvec::Vector)

    abstol, reltol = get_option(m, :disc_abs_width_tol), get_option(m, :disc_rel_width_tol)

    lb_local, ub_local = partvec[partidx], partvec[partidx+1]
    ub_touch, lb_touch = true, true
    lb_new, ub_new = max(point - radius, lb_local), min(point + radius, ub_local)

    if ub_new < ub_local && !isapprox(ub_new, ub_local; atol=abstol) && abs(ub_new-ub_local)/(1e-8+abs(ub_local)) > reltol # Insert new UB-based partition
        insert!(partvec, partidx+1, ub_new)
        ub_touch = false
    end

    if lb_new > lb_local && !isapprox(lb_new, lb_local; atol=abstol) && abs(lb_new-lb_local)/(1e-8+abs(lb_local)) > reltol # Insert new LB-based partition
        insert!(partvec, partidx+1, lb_new)
        lb_touch = false
    end

    if (ub_touch && lb_touch) || check_solution_history(m, var)
        distvec = [(j, partvec[j+1]-partvec[j]) for j in 1:length(partvec)-1]
        sort!(distvec, by=x->x[2])
        point_orig = point
        pos = distvec[end][1]
        lb_local = partvec[pos]
        ub_local = partvec[pos+1]
        isapprox(lb_local, ub_local;atol=get_option(m, :tol)) && return
        chunk = (ub_local - lb_local) / get_option(m, :disc_divert_chunks)
        point = lb_local + (ub_local - lb_local) / get_option(m, :disc_divert_chunks)
        for i in 2:get_option(m, :disc_divert_chunks)
            insert!(partvec, pos+1, lb_local + chunk * (get_option(m, :disc_divert_chunks)-(i-1)))
        end
        (get_option(m, :loglevel) > 199) && println("[DEBUG] !D! VAR$(var): SOL=$(round(point_orig; digits=4))=>$(point) |$(round(lb_local; digits=4)) | $(get_option(m, :disc_divert_chunks)) SEGMENTS | $(round(ub_local; digits=4))|")
    else
        (get_option(m, :loglevel) > 199) && println("[DEBUG] VAR$(var): SOL=$(round(point; digits=4)) RADIUS=$(radius), PARTITIONS=$(length(partvec)-1) |$(round(lb_local; digits=4)) |$(round(lb_new; digits=6)) <- * -> $(round(ub_new; digits=6))| $(round(ub_local; digits=4))|")
    end

    return
end

function add_uniform_partition(m::Optimizer; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization

    for i in m.disc_vars  # Only construct when discretized
        lb_local = discretization[i][1]
        ub_local = discretization[i][end]
        distance = ub_local - lb_local
        chunk = distance / ((m.logs[:n_iter]+1)*get_option(m, :disc_uniform_rate))
        discretization[i] = [lb_local+chunk*(j-1) for j in 1:(m.logs[:n_iter]+1)*get_option(m, :disc_uniform_rate)]
        push!(discretization[i], ub_local)   # Safety Scheme
        (get_option(m, :loglevel) > 199) && println("[DEBUG] VAR$(i): RATE=$(get_option(m, :disc_uniform_rate)), PARTITIONS=$(length(discretization[i]))  |$(round(lb_local; digits=4)) | $(get_option(m, :disc_uniform_rate)*(1+m.logs[:n_iter])) SEGMENTS | $(round(ub_local; digits=4))|")
    end

    return discretization
end

function update_disc_ratio(m::Optimizer, presolve=false)

    m.logs[:n_iter] > 2 && return get_option(m, :disc_ratio) # Stop branching after the second iterations

    ratio_pool = [8:2:20;]  # Built-in try range
    convertor = Dict(MOI.MAX_SENSE => :<, MOI.MIN_SENSE => :>)
    revconvertor = Dict(MOI.MAX_SENSE => :>, MOI.MIN_SENSE => :<)

    incumb_ratio = ratio_pool[1]
    is_min_sense(m) ? incumb_res = -Inf : incumb_res = Inf
    res_collector = Float64[]

    for r in ratio_pool
        st = time()
        if !isempty(m.best_sol)
            branch_disc = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=r, use_solution=m.best_sol)
        else
            branch_disc = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=r)
        end
        create_bounding_mip(m, use_disc=branch_disc)
        res = disc_branch_solve(m)
        push!(res_collector, res)
        if eval(convertor[m.sense_orig])(res, incumb_res) # && abs(abs(collector[end]-res)/collector[end]) > 1e-1  # %1 of difference
            incumb_res = res
            incumb_ratio = r
        end
        et = time() - st
        if et > 300  # 5 minutes limitation
            println("Expensive disc branching pass... Fixed at 8")
            return 8
        end
        get_option(m, :loglevel) > 0 && println("BRANCH RATIO = $(r), METRIC = $(res) || TIME = $(time()-st)")
    end

    if std(res_collector) >= 1e-2    # Detect if all solution are similar to each other
        get_option(m, :loglevel) > 0 && println("RATIO BRANCHING OFF due to solution variance test passed.")
        set_option(m, :disc_ratio_branch, false) # If an incumbent ratio is selected, then stop the branching scheme
    end

    if !isempty(m.best_sol)
        m.discretization = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=incumb_ratio, use_solution=m.best_sol)
    else
        m.discretization = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=incumb_ratio)
    end

    get_option(m, :loglevel) > 0 && println("INCUMB_RATIO = $(incumb_ratio)")

    return incumb_ratio
end

function disc_branch_solve(m::Optimizer)

    # ================= Solve Start ================ #
    update_mip_time_limit(m)
    start_bounding_solve = time()
    status = solve(m.model_mip, suppress_warnings=true)
    cputime_branch_bounding_solve = time() - start_bounding_solve
    m.logs[:total_time] += cputime_branch_bounding_solve
    m.logs[:time_left] = max(0.0, get_option(m, :timeout) - m.logs[:total_time])
    # ================= Solve End ================ #

    if status in [:Optimal, :Suboptimal, :UserLimit]
        return m.model_mip.objBound
    else
        @warn "  Warning: Unexpected solving condition $(status) during disc branching."
    end

    # Safety scheme
    if is_min_sense(m)
        return -Inf
    else
        return Inf
    end
end
