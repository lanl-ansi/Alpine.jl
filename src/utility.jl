"""
    update_rel_gap(m::PODNonlinearModel)

Update POD model relative & absolute optimality gap.

The relative gap calculation is

```math
    \\textbf{Gap} = \\frac{|UB-LB|}{ϵ+|UB|}
```

The absolute gap calculation is
```
    |UB-LB|
```
"""
function update_opt_gap(m::PODNonlinearModel)

    m.best_rel_gap = abs(m.best_obj - m.best_bound)/(m.tol+abs(m.best_obj))
    # absoluate or anyother bound calculation shows here...

    return
end

"""
    update_var_bounds(m::PODNonlinearModel, discretization::Dict; len::Float64=length(keys(discretization)))

This function take in a dictionary-based discretization information and convert them into two bounds vectors (l_var, u_var) by picking the smallest and largest numbers. User can specify a certain length that may contains variables that is out of the scope of discretization.

Output::

    l_var::Vector{Float64}, u_var::Vector{Float64}
"""
function update_var_bounds(discretization; kwargs...)

    options = Dict(kwargs)

    haskey(options, :len) ? len = options[:len] : len = length(keys(discretization))

    l_var = fill(-Inf, len)
    u_var = fill(Inf, len)

    for var_idx in keys(discretization)
        l_var[var_idx] = discretization[var_idx][1]
        u_var[var_idx] = discretization[var_idx][end]
    end

    return l_var, u_var
end

"""
    discretization_to_bounds(d::Dict, l::Int)

Same as [`update_var_bounds`](@ref)
"""
discretization_to_bounds(d::Dict, l::Int) = update_var_bounds(d, len=l)

"""
    initialize_discretization(m::PODNonlinearModel)

This function initialize the dynamic discretization used for any bounding models. By default, it takes (.l_var_orig, .u_var_orig) as the base information. User is allowed to use alternative bounds for initializing the discretization dictionary.
The output is a dictionary with MathProgBase variable indices keys attached to the :PODNonlinearModel.discretization.
"""
function initialize_discretization(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    for var in 1:(m.num_var_orig+m.num_var_lifted_mip)
        lb = m.l_var_tight[var]
        ub = m.u_var_tight[var]
        m.discretization[var] = [lb, ub]
    end

    return
end

"""

    to_discretization(m::PODNonlinearModel, lbs::Vector{Float64}, ubs::Vector{Float64})

Utility functions to convert bounds vectors to Dictionary based structures that is more suitable for
partition operations.

"""
function to_discretization(m::PODNonlinearModel, lbs::Vector{Float64}, ubs::Vector{Float64}; kwargs...)

    options = Dict(kwargs)

    var_discretization = Dict()
    for var in 1:m.num_var_orig
        lb = lbs[var]
        ub = ubs[var]
        var_discretization[var] = [lb, ub]
    end

    for var in (1+m.num_var_orig):(m.num_var_orig+m.num_var_lifted_mip)
        lb = -Inf
        ub = Inf
        var_discretization[var] = [lb, ub]
    end

    return var_discretization
end

"""
    flatten_discretization(discretization::Dict)

Utility functions to eliminate all partition on discretizing variable and keep the loose bounds.

"""
function flatten_discretization(discretization::Dict; kwargs...)

    flatten_discretization = Dict()
    for var in keys(discretization)
        flatten_discretization[var] = [discretization[var][1],discretization[var][end]]
    end

    return flatten_discretization
end

"""
    add_discretization(m::PODNonlinearModel; use_discretization::Dict, use_solution::Vector)

Basic built-in method used to add a new partition on feasible domains of discretizing variables.
This method make modification in .discretization

Consider original partition [0, 3, 7, 9], where LB/any solution is 4.
Use ^ as the new partition, "|" as the original partition

A case when discretize ratio = 4
| -------- | - ^ -- * -- ^ ---- | -------- |
0          3  3.5   4   4.5     7          9

A special case when discretize ratio = 2
| -------- | ---- * ---- ^ ---- | -------- |
0          3      4      5      7          9

There are two options for this function,

    * `use_discretization(default=m.discretization)`:: to regulate which is the base to add new partitions on
    * `use_solution(default=m.best_bound_sol)`:: to regulate which solution to use when adding new partitions on

TODO: also need to document the speical diverted cases when new partition touches both sides

This function belongs to the hackable group, which means it can be replaced by the user to change the behvaior of the solver.
"""
function add_discretization(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization
    haskey(options, :use_solution) ? point_vec = options[:use_solution] : point_vec = m.best_bound_sol

    # ? Perform discretization base on type of nonlinear terms
    for i in 1:m.num_var_orig
        point = point_vec[i]
        # @show i, point, discretization[i]
        @assert point >= discretization[i][1] - m.tol       # Solution validation
        @assert point <= discretization[i][end] + m.tol
        # Safety Scheme
        (abs(point - discretization[i][1]) <= m.tol) && (point = discretization[i][1])
        (abs(point - discretization[i][end]) <= m.tol) && (point = discretization[i][end])
        if i in m.var_discretization_mip  # Only construct when discretized
            for j in 1:length(discretization[i])
                if point >= discretization[i][j] && point <= discretization[i][j+1]  # Locating the right location
                    @assert j < length(m.discretization[i])
                    lb_local = discretization[i][j]
                    ub_local = discretization[i][j+1]
                    distance = ub_local - lb_local
                    radius = distance / m.discretization_ratio
                    lb_new = max(point - radius, lb_local)
                    ub_new = min(point + radius, ub_local)
                    ub_touch = true
                    lb_touch = true
                    if ub_new < ub_local && !isapprox(ub_new, ub_local; atol=1e-3)  # Insert new UB-based partition
                        insert!(discretization[i], j+1, ub_new)
                        ub_touch = false
                    end
                    if lb_new > lb_local && !isapprox(lb_new, lb_local; atol=1e-3)  # Insert new LB-based partition
                        insert!(discretization[i], j+1, lb_new)
                        lb_touch = false
                    end
                    # if ub_touch && lb_touch
                    #     distance = -1.0
                    #     pos = -1
                    #     for j in 2:length(discretization[i])  # it is made sure there should be at least two partitions
                    #         if (discretization[i][j] - discretization[i][j-1]) > distance
                    #             lb_local = discretization[i][j-1]
                    #             ub_local = discretization[i][j]
                    #             distance = ub_local - lb_local
                    #             point = lb_local + (ub_local - lb_local) / 2   # reset point
                    #             pos = j
                    #         end
                    #     end
                    #     rdius = distance / m.discretization_ratio
                    #     lb_new = max(point - radius, lb_local)
                    #     ub_new = min(point + radius, ub_local)
                    #     if ub_new < ub_local && !isapprox(ub_new, ub_local; atol=m.tol)  # Insert new UB-based partition
                    #         insert!(discretization[i], pos, ub_new)
                    #     end
                    #     if lb_new > lb_local && !isapprox(lb_new, lb_local; atol=m.tol)  # Insert new LB-based partition
                    #         insert!(discretization[i], pos, lb_new)
                    #     end
                    #     m.log_level > 99 && println("[DEBUG] VAR$(i): !diverted! : SOL=$(round(point,4)) RATIO=$(m.discretization_ratio), PARTITIONS=$(length(discretization[i])-1)  |$(round(lb_local,4)) |$(round(lb_new,6)) <- * -> $(round(ub_new,6))| $(round(ub_local,4))|")
                    # else
                    #     m.log_level > 99 && println("[DEBUG] VAR$(i): SOL=$(round(point,4)) RATIO=$(m.discretization_ratio), PARTITIONS=$(length(discretization[i])-1)  |$(round(lb_local,4)) |$(round(lb_new,6)) <- * -> $(round(ub_new,6))| $(round(ub_local,4))|")
                    # end
                    break
                end
            end
        end
    end

    return discretization
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_mip_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    for i in 1:length(m.mip_solver.options)
        if fetch_timeleft_symbol(m) in collect(m.mip_solver.options[i])
            deleteat!(m.mip_solver.options, i)
            break
        end
    end

    if m.timeout != Inf
        push!(m.mip_solver.options, (fetch_timeleft_symbol(m), timelimit))
    end

    return
end

"""

    fetch_timeleft_symbol(m::PODNonlinearModel)

An utility function used to recognize differnt sub-solvers return the timelimit setup keywords.
"""
function fetch_timeleft_symbol(m::PODNonlinearModel)
    if string(m.mip_solver)[1:5] == "CPLEX"
        return :CPX_PARAM_TILIM
    elseif string(m.mip_solver)[1:6] == "Gurobi"
        return :TimeLimit
    elseif string(m.mip_solver)[1:3] == "Cbc"
        return :seconds
    else found == nothing
        error("Needs support for this MIP solver")
    end
    return
end

"""
    fetch_boundstop_symbol(m::PODNonlinearModel)

An utility function used to recongize different sub-solvers and return the bound stop option key words
"""
function update_boundstop_options(m::PODNonlinearModel)

    # Calculation of the bound
    if m.sense_orig == :Min
        stopbound = (1-m.rel_gap+m.tol) * m.best_obj
    elseif m.sense_orig == :Max
        stopbound = (1+m.rel_gap-m.tol) * m.best_obj
    end

    for i in 1:length(m.mip_solver.options)
        if m.mip_solver.options[i][1] == :BestBdStop
            deleteat!(m.mip_solver.options, i)
            if string(m.mip_solver)[1:6] == "Gurobi"
                push!(m.mip_solver.options, (:BestBdStop, stopbound))
            else
                return
            end
        end
    end

    if string(m.mip_solver)[1:6] == "Gurobi"
        push!(m.mip_solver.options, (:BestBdStop, stopbound))
    else
        return
    end

    return
end

"""

    fix_domains(m::PODNonlinearModel)

This function is used to fix variables to certain domains during the local solve process in the [`global_solve`](@ref).
More specifically, it is used in [`local_solve`](@ref) to fix binary and integer variables to lower bound solutions
and discretizing varibles to the active domain according to lower bound solution.
"""
function fix_domains(m::PODNonlinearModel; kwargs...)

    l_var = copy(m.l_var_orig)
    u_var = copy(m.u_var_orig)
    for i in 1:m.num_var_orig
        if i in m.var_discretization_mip
            point = m.sol_incumb_lb[i]
            for j in 1:length(m.discretization[i])
                if point >= (m.discretization[i][j] - m.tol) && (point <= m.discretization[i][j+1] + m.tol)
                    @assert j < length(m.discretization[i])
                    l_var[i] = m.discretization[i][j]
                    u_var[i] = m.discretization[i][j+1]
                    break
                end
            end
        elseif m.var_type_orig[i] == :Bin || m.var_type_orig[i] == :Int
            l_var[i] = round(m.sol_incumb_lb[i])
            u_var[i] = round(m.sol_incumb_lb[i])
        end
    end

    return l_var, u_var
end

"""
    convexification_exam(m::PODNonlinearModel)
"""
function convexification_exam(m::PODNonlinearModel)

    # Other more advanced convexification check goes here
    for term in keys(m.nonlinear_info)
        if !m.nonlinear_info[term][:convexified]
            warn("Detected terms that is not convexified $(term[:lifted_constr_ref]), bounding model solver may report a error due to this")
            return
        end
    end

    return
end

"""
    pick_vars_discretization(m::PODNonlinearModel)

This function helps pick the variables for discretization. The method chosen depends on user-inputs.
In case when `indices::Int` is provided, the method is chosen as built-in method. Currently,
there exist two built-in method:

    * `max-cover(m.discretization_var_pick_algo=0, default)`: pick all variables involved in the non-linear term for discretization
    * `min-vertex-cover(m.discretization_var_pick_algo=1)`: pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified

For advance usage, `m.discretization_var_pick_algo` allows `::Function` inputs. User is required to perform flexible methods in choosing the non-linear variable.
For more information, read more details at [Hacking Solver](@ref).

"""
function pick_vars_discretization(m::PODNonlinearModel)

    if isa(m.discretization_var_pick_algo, Function)
        (m.log_level > 0) && println("using method $(m.discretization_var_pick_algo) for picking discretization variable...")
        eval(m.discretization_var_pick_algo)(m)
        (length(m.var_discretization_mip) == 0 && length(m.nonlinear_info) > 0) && error("[USER FUNCTION] must select at least one variable to perform discretization for convexificiation purpose")
    elseif isa(m.discretization_var_pick_algo, Int) || isa(m.discretization_var_pick_algo, String)
        if m.discretization_var_pick_algo == 0 || m.discretization_var_pick_algo == "max_cover"
            max_cover(m)
        elseif m.discretization_var_pick_algo == 1 || m.discretization_var_pick_algo == "min_vertex_cover"
            min_vertex_cover(m)
        else
            error("Unsupported default indicator for picking variables for discretization")
        end
    else
        error("Input for parameter :discretization_var_pick_algo is illegal. Should be either a Int for default methods indexes or functional inputs.")
    end

    return
end

"""
    mccormick(m::JuMP.Model, xy, x, y, x_l, x_u, y_l, y_u)

Generic function to add a McCormick convex envelop, where `xy=x*y` and `x_l, x_u, y_l, y_u` are variable bounds.
"""
function mccormick(m,xy,x,y,xˡ,xᵘ,yˡ,yᵘ)
    @constraint(m, xy >= xˡ*y + yˡ*x - xˡ*yˡ)
    @constraint(m, xy >= xᵘ*y + yᵘ*x - xᵘ*yᵘ)
    @constraint(m, xy <= xˡ*y + yᵘ*x - xˡ*yᵘ)
    @constraint(m, xy <= xᵘ*y + yˡ*x - xᵘ*yˡ)
    return
end

function mccormick_bin(m,xy,x,y)

    @constraint(m, xy <= x)
    @constraint(m, xy <= y)
    @constraint(m, xy >= x+y-1)

    return
end

function mccormick_monomial(m,xy,x,xˡ,xᵘ)
    @constraint(m, xy >= x^2)
    @constraint(m, xy <= (xˡ+xᵘ)*x - (xˡ*xᵘ))
    return
end

function tightmccormick_monomial(m,x_p,x,xz,xˡ,xᵘ,z,p,lazy,quad) # if p=2, tightened_lazycuts = tightmccormick_quad
    if lazy == 1
        function GetLazyCuts_quad(cb)
            TOL1 = 1e-6
            if (getvalue(x)^p > (getvalue(x_p) + TOL1))
                a = p*getvalue(x)^(p-1)
                b = (1-p)*getvalue(x)^p
                @lazyconstraint(cb, a*x + b <= x_p)
            end
        end
        addlazycallback(m, GetLazyCuts_quad)
    elseif p == 2 && quad == 1
        @constraint(m, x_p >= x^2)
    else
        x0_vec = sort(union(xˡ, xᵘ))
        for x0 in x0_vec
            @constraint(m, x_p >= (1-p)*(x0)^p + p*(x0)^(p-1)*x)
        end
    end

    A = ((xᵘ).^p-(xˡ).^p)./(xᵘ-xˡ)
    @constraint(m, x_p .<= A'*xz - (A.*xˡ)'*z + ((xˡ).^p)'*z)

    return
end

"""

    min_vertex_cover(m:PODNonlinearModel)

A built-in method for selecting variables for discretization.

"""
function min_vertex_cover(m::PODNonlinearModel)

    # Collect the information for arcs and nodes
    nodes = Set()
    arcs = Set()
    for pair in keys(m.nonlinear_info)
        arc = []
        if length(pair) > 2
            warn("min_vertex_cover discretizing variable selection method only support bi-linear problems, enfocing thie method may produce mistakes...")
        end
        for i in pair
            @assert isa(i.args[2], Int)
            push!(nodes, i.args[2])
            push!(arc, i.args[2])
        end
        push!(arcs, arc)
    end
    nodes = collect(nodes)
    arcs = collect(arcs)

    # Set up minimum vertex cover problem
    minvertex = Model(solver=m.mip_solver)
    @variable(minvertex, x[nodes], Bin)
    for arc in arcs
        @constraint(minvertex, x[arc[1]] + x[arc[2]] >= 1)
    end
    @objective(minvertex, Min, sum(x))
    status = solve(minvertex, suppress_warnings=true)

    xVal = getvalue(x)

    # Getting required information
    m.num_var_discretization_mip = Int(sum(xVal))
    m.var_discretization_mip = [i for i in nodes if xVal[i] > 1e-5]

end

"""

    max_cover(m:PODNonlinearModel)

A built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.

"""
function max_cover(m::PODNonlinearModel; kwargs...)

    nodes = Set()
    for pair in keys(m.nonlinear_info)
        if length(pair) > 2
            warn("max-cover discretizing variable selection method only support bi-linear problems. enforcing may produce mistakes...")
        end
        for i in pair
            @assert isa(i.args[2], Int)
            push!(nodes, i.args[2])
        end
    end
    nodes = collect(nodes)
    m.num_var_discretization_mip = length(nodes)
    m.var_discretization_mip = nodes

    return
end
