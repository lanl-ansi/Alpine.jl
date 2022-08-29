"""
   update_opt_gap(m::Optimizer)

Updates Alpine model's relative & absolute optimality gaps.

The relative gap calculation is

```math
\\textbf{Gap} = \\frac{|UB-LB|}{ϵ+|UB|}
```

The absolute gap calculation is
```
|UB-LB|
```
"""
function update_opt_gap(m::Optimizer)
    tol = Alp.get_option(m, :tol)

    if m.best_obj in [Inf, -Inf]
        m.best_rel_gap = Inf
        return
    else
        p = convert(Int, round(abs(log(10, Alp.get_option(m, :rel_gap)))))
        n = round(abs(m.best_obj - m.best_bound); digits = p)
        # dn = round(abs(1e-12+abs(m.best_obj)); digits=p)
        if isapprox(n, 0.0; atol = tol) && isapprox(m.best_obj, 0.0; atol = tol)
            m.best_rel_gap = 0.0
            return
        end

        if Alp.is_min_sense(m)
            m.best_rel_gap = Alp.eval_opt_gap(m, m.best_bound, m.best_obj)
        elseif Alp.is_max_sense(m)
            m.best_rel_gap = Alp.eval_opt_gap(m, m.best_obj, m.best_bound)
        end
    end

    m.best_abs_gap = abs(m.best_obj - m.best_bound)

    return
end

function eval_opt_gap(m::Optimizer, lower_bound::Number, upper_bound::Number)
    tol = Alp.get_option(m, :tol)

    if isapprox(lower_bound, upper_bound, atol = tol)
        m.best_rel_gap = 0.0

    elseif (lower_bound - upper_bound) > 1E-5
        error("Lower bound cannot cross the upper bound in optimality gap evaluation")

    elseif isinf(lower_bound) || isinf(upper_bound)
        error("Infinite bounds detected in optimality gap evalutation")

    else
        if isapprox(upper_bound, 0.0; atol = tol) # zero upper bound case
            eps = 1 # shift factor
            m.best_rel_gap =
                abs((upper_bound + eps) - (lower_bound + eps)) /
                (tol + abs(upper_bound) + eps)
        else
            m.best_rel_gap = abs(upper_bound - lower_bound) / (tol + abs(upper_bound))
        end
    end

    return m.best_rel_gap
end

function measure_relaxed_deviation(m::Optimizer; sol = nothing)
    sol = something(sol, m.best_bound_sol)

    isempty(sol) && return

    dev = []
    for k in keys(m.nonconvex_terms)
        y_idx = m.nonconvex_terms[k][:y_idx]
        y_hat = sol[y_idx]
        y_val = m.nonconvex_terms[k][:evaluator](m.nonconvex_terms[k], sol)
        push!(
            dev,
            (y_idx, abs(y_hat - y_val), y_hat, y_val, m.nonconvex_terms[k][:var_idxs]),
        )
    end

    sort!(dev, by = x -> x[1])

    for i in dev
        Alp.get_option(m, :log_level) > 199 && println(
            "Y-VAR$(i[1]): DIST=$(i[2]) || Y-hat = $(i[3]), Y-val = $(i[4]) || COMP $(i[5])",
        )
    end

    return
end

"""
   discretization_to_bounds(d::Dict, l::Int)

Same as [`update_var_bounds`](@ref)
"""
discretization_to_bounds(d::Dict, l::Int) = Alp.update_var_bounds(d, len = l)

"""
Update the data structure with feasible solution and its associated objective (if better)
"""
function update_incumb_objective(m::Optimizer, objval::Float64, sol::Vector)
    convertor = Dict(MOI.MAX_SENSE => :>, MOI.MIN_SENSE => :<)
    push!(m.logs[:obj], objval)
    if eval(convertor[m.sense_orig])(objval, m.best_obj) #&& !eval(convertor[m.sense_orig])(objval, m.best_bound)
        m.best_obj = objval
        m.best_sol = sol
        m.detected_feasible_solution = true
    end

    return
end

"""
   variable_values(m::JuMP.Model)

Prints values of all the variables after optimizing the original JuMP model.
"""
function variable_values(m::JuMP.Model)
    for var in all_variables(m)
        if isapprox(abs(JuMP.value(var)), 0, atol = 1E-6)
            println("$var = 0.0")
        else
            println("$var = $(round(JuMP.value(var), digits = 5))")
        end
    end

    return
end

"""
Mention the variable type of the lifted term.
This function is with limited functionality

"""
function resolve_lifted_var_type(var_types::Vector{Symbol}, operator::Symbol)
    if operator == :+
        detector = [i in [:Bin, :Int] ? true : false for i in var_types]
        if length(detector) == 1 && detector[1] # Special case
            if :Bin in var_types
                return :Bin
            else
                return :Int
            end
        end
        prod(detector) && return :Int
        # o/w continous variables

    elseif operator == :*
        detector = [i == :Bin ? true : false for i in var_types]
        prod(detector) && return :Bin
        detector = [i in [:Bin, :Int] ? true : false for i in var_types]
        prod(detector) && return :Int
        # o/w continous variables

    end

    return :Cont
end

"""
   check_solution_history(m::Optimizer, ind::Int)

Check if the solution is always the same within the last disc_consecutive_forbid iterations. Return `true` if solution has stalled.
"""
function check_solution_history(m::Optimizer, ind::Int)
    Alp.get_option(m, :disc_consecutive_forbid) == 0 && return false
    (m.logs[:n_iter] < Alp.get_option(m, :disc_consecutive_forbid)) && return false

    sol_val = m.bound_sol_history[mod(
        m.logs[:n_iter] - 1,
        Alp.get_option(m, :disc_consecutive_forbid),
    )+1][ind]
    for i in 1:(Alp.get_option(m, :disc_consecutive_forbid)-1)
        search_pos =
            mod(m.logs[:n_iter] - 1 - i, Alp.get_option(m, :disc_consecutive_forbid)) + 1
        !isapprox(
            sol_val,
            m.bound_sol_history[search_pos][ind];
            atol = Alp.get_option(m, :disc_rel_width_tol),
        ) && return false
    end

    Alp.get_option(m, :log_level) > 99 &&
        println("Consecutive bounding solution on VAR$(ind) obtained. Diverting...")

    return true
end

"""
   fix_domains(m::Optimizer)

This function is used to fix variables to certain domains during the local solve process in the [`global_solve`](@ref).
More specifically, it is used in [`local_solve`](@ref) to fix binary and integer variables to lower bound solutions
and discretizing variables to the active domain according to lower bound solution.
"""
function fix_domains(m::Optimizer; discrete_sol = nothing, use_orig = false)
    discrete_sol !== nothing && @assert length(discrete_sol) >= m.num_var_orig

    l_var = [m.l_var_tight[i] for i in 1:m.num_var_orig]
    u_var = [m.u_var_tight[i] for i in 1:m.num_var_orig]

    for i in 1:m.num_var_orig
        if i in m.disc_vars && m.var_type[i] == :Cont
            point = if discrete_sol === nothing
                m.best_bound_sol[i]
            else
                discrete_sol[i]
            end
            PCnt = length(m.discretization[i]) - 1
            for j in 1:PCnt
                if point >= (m.discretization[i][j] - Alp.get_option(m, :tol)) &&
                   (point <= m.discretization[i][j+1] + Alp.get_option(m, :tol))
                    @assert j < length(m.discretization[i])
                    use_orig ? l_var[i] = m.discretization[i][1] :
                    l_var[i] = m.discretization[i][j]
                    use_orig ? u_var[i] = m.discretization[i][end] :
                    u_var[i] = m.discretization[i][j+1]
                    break
                end
            end
        elseif m.var_type[i] == :Bin || m.var_type[i] == :Int
            if discrete_sol === nothing
                l_var[i] = round(m.best_bound_sol[i])
                u_var[i] = round(m.best_bound_sol[i])
            else
                l_var[i] = round(discrete_sol[i])
                u_var[i] = round(discrete_sol[i])
            end
        end
    end

    return l_var, u_var
end

"""
   is_fully_convexified(m::Optimizer)
"""
function is_fully_convexified(m::Optimizer)

    # Other more advanced convexification check goes here
    for term in keys(m.nonconvex_terms)
        if !m.nonconvex_terms[term][:convexified]
            @warn "  Warning: Detected terms that is not convexified $(term[:lifted_constr_ref]), bounding model solver may report a error due to this"
            return
        else
            m.nonconvex_terms[term][:convexified] = false    # Reset status for next iteration
        end
    end

    return
end

"""
Collect active partition idx
Need to be careful with the diverted point
"""
function get_active_partition_idx(
    discretization::Dict,
    val::Float64,
    idx::Int;
    tol = 1e-6,
)
    for j in 1:length(discretization[idx])-1
        if val > discretization[idx][j] - tol && val < discretization[idx][j+1] + tol
            return j
        end
    end

    @warn "  Warning: Activate parition not found [VAR$(idx)]. Returning default partition 1."
    return 1
end

"""
   get_candidate_disc_vars(m:Optimizer)

   A built-in method for selecting variables for discretization. This function selects all candidate variables in the nonlinear terms for discretization,
   if under a threshold value of the number of nonlinear terms.
"""
function get_candidate_disc_vars(m::Optimizer; getoutput = false)

    # Pick variables that is bound width more than tolerance length
    if getoutput
        return [i for i in m.candidate_disc_vars]
    else
        m.disc_vars = [i for i in m.candidate_disc_vars]
        m.num_var_disc_mip = length(m.disc_vars)
    end

    return
end

function initialize_solution_pool(m::Optimizer, cnt::Int)
    s = Dict()

    s[:cnt] = cnt

    # Column dimension changing variable
    s[:len] = m.num_var_orig + m.num_var_linear_mip + m.num_var_nonlinear_mip

    # !! Be careful with the :vars when utilizing the dynamic discretization variable selection !!
    s[:vars] = [i for i in m.candidate_disc_vars if length(m.discretization[i]) > 2]

    s[:sol] = Vector{Vector}(undef, cnt)                   # Solution value
    s[:obj] = Vector{Float64}(undef, cnt)                  # Objecitve value
    s[:disc] = Vector{Dict}(undef, cnt)                    # Discretization
    s[:stat] = [:Alive for i in 1:cnt]              # Solution status
    s[:iter] = [m.logs[:n_iter] for i in 1:cnt]     # Iteration collected
    s[:ubstart] = [false for i in 1:cnt]            # Solution used for ub multistart

    return s
end

"""
Reconsideration required
"""
function ncvar_collect_arcs(m::Optimizer, nodes::Vector)
    arcs = Set()

    for k in keys(m.nonconvex_terms)
        if m.nonconvex_terms[k][:nonlinear_type] == :BILINEAR
            arc = [i.args[2] for i in k]
            length(arc) == 2 && push!(arcs, sort(arc))

        elseif m.nonconvex_terms[k][:nonlinear_type] == :MONOMIAL
            @assert isa(m.nonconvex_terms[k][:var_idxs][1], Int)
            varidx = m.nonconvex_terms[k][:var_idxs][1]
            push!(arcs, [varidx; varidx])

        elseif m.nonconvex_terms[k][:nonlinear_type] == :MULTILINEAR
            varidxs = m.nonconvex_terms[k][:var_idxs]
            for i in 1:length(varidxs)
                for j in 1:length(varidxs)
                    if i != j
                        push!(arcs, sort([varidxs[i]; varidxs[j]]))
                    end
                end
            end
            if length(varidxs) == 1
                push!(arcs, sort([varidxs[1]; varidxs[1]]))
            end

            # elseif m.nonconvex_terms[k][:nonlinear_type] == :INTLIN
            #    var_idxs = copy(m.nonconvex_terms[k][:var_idxs])
            #    push!(arcs, sort(var_idxs))
            # elseif m.nonconvex_terms[k][:nonlinear_type] == :INTPROD
            #    var_idxs = m.nonconvex_terms[k][:var_idxs]
            #    for i in 1:length(var_idxs)
            #       for j in 1:length(var_idxs)
            #          i != j && push!(arcs, sort([var_idxs[i]; var_idxs[j]]))
            #       end
            #    end
            # elseif m.nonconvex_terms[k][:nonlinear_type] in [:cos, :sin]
            #    @assert length(m.nonconvex_terms[k][:var_idxs]) == 1
            #    var_idx = m.nonconvex_terms[k][:var_idxs][1]
            #    push!(arcs, [var_idx; var_idx])

        elseif m.nonconvex_terms[k][:nonlinear_type] in [:BININT, :BINLIN, :BINPROD]
            continue
        else
            error(
                "[EXCEPTION] Unexpected nonlinear term when building interaction graphs for min. vertex cover.",
            )
        end
    end

    return arcs
end

"""
"""
function build_discvar_graph(m::Optimizer)

    # Collect the information of nonlinear terms in terms of arcs and nodes
    nodes = Alp.get_candidate_disc_vars(m, getoutput = true)
    arcs = Alp.ncvar_collect_arcs(m, nodes)

    # Collect integer variables
    for i in 1:m.num_var_orig
        if !(i in nodes) && m.var_type[i] == :Int
            push!(nodes, i)
            push!(arcs, [i; i])
        end
    end

    nodes = collect(nodes)
    arcs = collect(arcs)

    return nodes, arcs
end

"""
   min_vertex_cover(m::Optimizer)

`min_vertex_cover` chooses the variables based on the minimum vertex cover algorithm for the interaction graph of
nonlinear terms which are adaptively partitioned for global optimization. This option can be activated by
setting `disc_var_pick = 1`.
"""
function min_vertex_cover(m::Optimizer)
    nodes, arcs = Alp.build_discvar_graph(m)

    # Set up minimum vertex cover problem
    minvertex = Model(Alp.get_option(m, :mip_solver))
    MOI.set(minvertex, MOI.TimeLimitSec(), 60.0) # Time limit for min vertex cover formulation
    JuMP.@variable(minvertex, 0 <= x[nodes] <= 1, Bin)
    JuMP.@constraint(minvertex, [a in arcs], x[a[1]] + x[a[2]] >= 1)
    JuMP.@objective(minvertex, Min, sum(x))

    JuMP.optimize!(minvertex)
    # status = MOI.get(minvertex, MOI.TerminationStatus())
    xVal = JuMP.value.(x)

    # Collecting required information
    m.num_var_disc_mip = Int(sum(xVal))
    m.disc_vars = [
        i for i in nodes if xVal[i] > Alp.get_option(m, :tol) &&
        abs(m.u_var_tight[i] - m.l_var_tight[i]) >= Alp.get_option(m, :tol)
    ]

    return
end

function weighted_min_vertex_cover(m::Optimizer, distance::Dict)

    # Collect the graph information
    nodes, arcs = Alp.build_discvar_graph(m)

    # A little bit redundancy before
    disvec = [distance[i] for i in keys(distance) if i in m.candidate_disc_vars]
    disvec = abs.(disvec[disvec.>0.0])
    isempty(disvec) ? heavy = 1.0 : heavy = 1 / minimum(disvec)
    weights = Dict()
    for i in m.candidate_disc_vars
        isapprox(distance[i], 0.0; atol = 1e-6) ? weights[i] = heavy :
        (weights[i] = (1 / distance[i]))
        (Alp.get_option(m, :log_level) > 100) &&
            println("VAR$(i) WEIGHT -> $(weights[i]) ||| DISTANCE -> $(distance[i])")
    end

    # Set up minimum vertex cover problem
    minvertex = Model(Alp.get_option(m, :mip_solver))
    MOI.set(minvertex, MOI.TimeLimitSec(), 60.0)  # Set a timer to avoid waste of time in proving optimality
    JuMP.@variable(minvertex, 0 <= x[nodes] <= 1, Bin)
    for arc in arcs
        JuMP.@constraint(minvertex, x[arc[1]] + x[arc[2]] >= 1)
    end
    JuMP.@objective(minvertex, Min, sum(weights[i] * x[i] for i in nodes))

    # Solve the minimum vertex cover
    JuMP.optimize!(minvertex)

    xVal = JuMP.value.(x)
    m.num_var_disc_mip = Int(sum(xVal))
    m.disc_vars = [
        i for i in nodes if xVal[i] > 0 &&
        abs(m.u_var_tight[i] - m.l_var_tight[i]) >= Alp.get_option(m, :tol)
    ]
    Alp.get_option(m, :log_level) >= 99 &&
        println("UPDATED DISC-VAR COUNT = $(length(m.disc_vars)) : $(m.disc_vars)")

    return
end

function round_sol(m::Optimizer, relaxed_sol)
    rounded_sol = copy(relaxed_sol)
    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Bin
            relaxed_sol[i] >= 0.5 ? rounded_sol[i] = 1 : rounded_sol[i] = 0
        elseif m.var_type_orig[i] == :Int
            rounded_sol[i] = round(relaxed_sol[i])
        else
            rounded_sol[i] = relaxed_sol[i]
        end
    end

    return rounded_sol
end

function _fetch_mip_solver_identifier(m::Optimizer; override="")
    (Alp.get_option(m, :mip_solver) === nothing) && return
    m.mip_solver_id = _get_solver_name(m, :mip_solver, override)
    return
end

function _fetch_nlp_solver_identifier(m::Optimizer; override="")
    (Alp.get_option(m, :nlp_solver) === nothing) && return
    m.nlp_solver_id = _get_solver_name(m, :nlp_solver, override)
    return
end

function _fetch_minlp_solver_identifier(m::Optimizer; override="")
    (Alp.get_option(m, :minlp_solver) === nothing) && return
    m.minlp_solver_id = _get_solver_name(m, :minlp_solver, override)
    return
end

function _get_solver_name(m::Optimizer, key::Symbol, override::String)
    solver_string = override
    if isempty(override)
        solver = MOI.instantiate(Alp.get_option(m, key))
        solver_string = try
            MOI.get(solver, MOI.SolverName())
        catch
            "$(typeof(solver))"
        end
    end
    for name in ("Ipopt", "Gurobi", "CPLEX", "Juniper")
        if occursin(name, solver_string)
            return name
        end
    end
    return solver_string
end

"""
    set_mip_time_limit(m::Optimizer)

An utility function used to dynamically regulate MILP solver time limits to fit Alpine solver's time limits.
"""
function set_mip_time_limit(m::Optimizer)
    time_limit = max(0.0, Alp.get_option(m, :time_limit) - m.logs[:total_time])
    return MOI.set(m.model_mip, MOI.TimeLimitSec(), time_limit)
end

"""
Follow the definition of terms to calculate the value of lifted terms
"""
function resolve_lifted_var_value(m::Optimizer, sol_vec::Array)
    @assert length(sol_vec) == m.num_var_orig
    sol_vec = [sol_vec; fill(NaN, m.num_var_linear_mip + m.num_var_nonlinear_mip)]

    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            lvar_idx = m.nonconvex_terms[k][:y_idx]
            sol_vec[lvar_idx] =
                m.nonconvex_terms[k][:evaluator](m.nonconvex_terms[k], sol_vec)
        elseif haskey(m.linear_terms, k)
            lvar_idx = m.linear_terms[k][:y_idx]
            sol_vec[lvar_idx] = m.linear_terms[k][:evaluator](m.linear_terms[k], sol_vec)
        else
            error("[RARE] Found homeless term key $(k) during bound resolution.")
        end
    end

    return sol_vec
end
