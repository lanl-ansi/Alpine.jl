"""
    update_rel_gap(m::AlpineNonlinearModel)

Update Alpine model relative & absolute optimality gap.

The relative gap calculation is

```math
    \\textbf{Gap} = \\frac{|UB-LB|}{ϵ+|UB|}
```

The absolute gap calculation is
```
    |UB-LB|
```
"""
function update_opt_gap(m::AlpineNonlinearModel)

    if m.best_obj in [Inf, -Inf]
        m.best_rel_gap = Inf
        return
    else
        p = convert(Int, round(abs(log(10,m.relgap))))
        n = round(abs(m.best_obj-m.best_bound); digits=p)
        dn = round(abs(1e-12+abs(m.best_obj)); digits=p)
        if isapprox(n, 0.0;atol=m.tol) && isapprox(m.best_obj,0.0;atol=m.tol)
            m.best_rel_gap = 0.0
            return
        end
        if m.gapref == :ub
            m.best_rel_gap = abs(m.best_obj - m.best_bound)/(m.tol+abs(m.best_obj))
        else
            m.best_rel_gap = abs(m.best_obj - m.best_bound)/(m.tol+abs(m.best_bound))
        end
    end

    m.best_abs_gap = abs(m.best_obj - m.best_bound)
    return
end

function measure_relaxed_deviation(m::AlpineNonlinearModel;sol=nothing)

    sol == nothing ? sol = m.best_bound_sol : sol = sol

    isempty(sol) && return

    dev = []
    for k in keys(m.nonconvex_terms)
        y_idx = m.nonconvex_terms[k][:y_idx]
        y_hat = sol[y_idx]
        y_val = m.nonconvex_terms[k][:evaluator](m.nonconvex_terms[k], sol)
        push!(dev, (y_idx, abs(y_hat-y_val), y_hat, y_val, m.nonconvex_terms[k][:var_idxs]))
    end

    sort!(dev, by=x->x[1])

    for i in dev
        m.loglevel > 199 && println("Y-VAR$(i[1]): DIST=$(i[2]) || Y-hat = $(i[3]), Y-val = $(i[4]) || COMP $(i[5])")
    end

    return
end

"""
    discretization_to_bounds(d::Dict, l::Int)

    Same as [`update_var_bounds`](@ref)
"""
discretization_to_bounds(d::Dict, l::Int) = update_var_bounds(d, len=l)

"""
    Update the data structure with feasible solution and its associated objective (if better)
"""
function update_incumb_objective(m::AlpineNonlinearModel, objval::Float64, sol::Vector)

    convertor = Dict(:Max=>:>, :Min=>:<)
    push!(m.logs[:obj], objval)
    if eval(convertor[m.sense_orig])(objval, m.best_obj) #&& !eval(convertor[m.sense_orig])(objval, m.best_bound)
        m.best_obj = objval
        m.best_sol = sol
        m.status[:feasible_solution] = :Detected
    end

    return
end

"""
    Utility function for debugging.
"""
function show_solution(m::JuMP.Model)
    for i in 1:length(m.colNames)
        println("$(m.colNames[i])=$(m.colVal[i])")
    end
    return
end


"""
    Tell what would be the variable type of a lifted term.
    This function is with limited functionality
    @docstring TODO
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
    @docstring
"""
function update_timeleft_symbol(options, keyword::Symbol, val::Float64; options_string_type=1)

    for i in 1:length(options)
        if options_string_type == 1
            if keyword in collect(options[i])
                deleteat!(options, i)
                break
            end
        elseif options_string_type == 2
            if keyword == split(options[i],"=")[1]
                deleteat!(options, i)
                break
            end
        end
    end

    if options_string_type == 1
        (val != Inf) && (push!(options, keyword => val))
    elseif options_string_type == 2
        (val != Inf) && (push!(options, "$(keyword)=$(val)"))
    end

    return options
end

"""
    fetch_boundstop_symbol(m::AlpineNonlinearModel)

An utility function used to recongize different sub-solvers and return the bound stop option key words
"""
function update_boundstop_options(m::AlpineNonlinearModel)

    if m.mip_solver_id == "Gurobi"
        # Calculation of the bound
        if m.sense_orig == :Min
            m.gapref == :ub ? stopbound=(1-m.relgap+m.tol)*abs(m.best_obj) : stopbound=(1-m.relgap+m.tol)*abs(m.best_bound)
        elseif m.sense_orig == :Max
            m.gapref == :ub ? stopbound=(1+m.relgap-m.tol)*abs(m.best_obj) : stopbound=(1+m.relgap-m.tol)*abs(m.best_bound)
        end

        for i in 1:length(m.mip_solver.options)
            if m.mip_solver.options[i][1] == :BestBdStop
                deleteat!(m.mip_solver.options, i)
                if m.mip_solver_id == "Gurobi"
                    push!(m.mip_solver.options, (:BestBdStop, stopbound))
                else
                    return
                end
            end
        end
    end

    return
end


"""
    check_solution_history(m::AlpineNonlinearModel, ind::Int)

Check if the solution is alwasy the same within the last disc_consecutive_forbid iterations. Return true if suolution in invariant.
"""
function check_solution_history(m::AlpineNonlinearModel, ind::Int)

    m.disc_consecutive_forbid == 0 && return false
    (m.logs[:n_iter] < m.disc_consecutive_forbid) && return false

    sol_val = m.bound_sol_history[mod(m.logs[:n_iter]-1, m.disc_consecutive_forbid)+1][ind]
    for i in 1:(m.disc_consecutive_forbid-1)
        search_pos = mod(m.logs[:n_iter]-1-i, m.disc_consecutive_forbid)+1
        !isapprox(sol_val, m.bound_sol_history[search_pos][ind]; atol=m.disc_rel_width_tol) && return false
    end

    m.loglevel > 99 && println("Consecutive bounding solution on VAR$(ind) obtained. Diverting...")
    return true
end

"""

    fix_domains(m::AlpineNonlinearModel)

This function is used to fix variables to certain domains during the local solve process in the [`global_solve`](@ref).
More specifically, it is used in [`local_solve`](@ref) to fix binary and integer variables to lower bound solutions
and discretizing varibles to the active domain according to lower bound solution.
"""
function fix_domains(m::AlpineNonlinearModel;discrete_sol=nothing, use_orig=false)

    discrete_sol != nothing && @assert length(discrete_sol) >= m.num_var_orig

    l_var = [m.l_var_tight[i] for i in 1:m.num_var_orig]
    u_var = [m.u_var_tight[i] for i in 1:m.num_var_orig]

    for i in 1:m.num_var_orig
        if i in m.disc_vars && m.var_type[i] == :Cont
            discrete_sol == nothing ? point = m.best_bound_sol[i] : point = discrete_sol[i]
            PCnt = length(m.discretization[i]) - 1
            for j in 1:PCnt
                if point >= (m.discretization[i][j] - m.tol) && (point <= m.discretization[i][j+1] + m.tol)
                    @assert j < length(m.discretization[i])
                    use_orig ? l_var[i] = m.discretization[i][1] : l_var[i] = m.discretization[i][j]
                    use_orig ? u_var[i] = m.discretization[i][end] : u_var[i] = m.discretization[i][j+1]
                    break
                end
            end
        elseif m.var_type[i] == :Bin || m.var_type[i] == :Int
            if discrete_sol == nothing
                l_var[i] = round(m.best_bound_sol[i])
                u_var[i] = round(m.best_bound_sol[i])
            else
                l_var[i] = discrete_sol[i]
                u_var[i] = discrete_sol[i]
            end
        end
    end

    return l_var, u_var
end

"""
    is_fully_convexified(m::AlpineNonlinearModel)
"""
function is_fully_convexified(m::AlpineNonlinearModel)

    # Other more advanced convexification check goes here
    for term in keys(m.nonconvex_terms)
        if !m.nonconvex_terms[term][:convexified]
            @warn "Detected terms that is not convexified $(term[:lifted_constr_ref]), bounding model solver may report a error due to this"
            return
        else
            m.nonconvex_terms[term][:convexified] = false    # Reset status for next iteration
        end
    end

    return
end

"""
    Collect LB solutions
    Don't test this function
"""
function collect_lb_pool(m::AlpineNonlinearModel)

    # Always stick to the structural .discretization for algorithm consideration info
    # If in need, the scheme need to be refreshed with customized discretization info

    if m.mip_solver_id != "Gurobi" || m.obj_structure == :convex || isempty([i for i in m.model_mip.colCat if i in [:Int, :Bin]])
        @warn "Skipping collecting solution pool procedure",
        return
    end

    # Construct a new solution pool for just this new iteration
    s = initialize_solution_pool(m, Gurobi.get_intattr(m.model_mip.internalModel.inner, "SolCount"))

    # Collect Solution and corresponding objective values
    for i in 1:s[:cnt]
        if m.mip_solver_id == "Gurobi"
            Gurobi.set_int_param!(m.model_mip.internalModel.inner, "SolutionNumber", i-1)
            s[:sol][i] = Gurobi.get_dblattrarray(m.model_mip.internalModel.inner, "Xn", 1, s[:len])
            s[:obj][i] = Gurobi.get_dblattr(m.model_mip.internalModel.inner, "PoolObjVal")
        elseif m.mip_solver_id == "CPLEX"
            error("No implementation for CPLEX")
        end
        s[:disc][i] = Dict(j=>get_active_partition_idx(m.discretization, s[:sol][i][j],j) for j in s[:vars])
    end

    merge_solution_pool(m, s)

    return
end

"""
    Merge collected solution pools
"""
function merge_solution_pool(m::AlpineNonlinearModel, s::Dict)

    # Always stick to the structural .discretization for algorithm consideration info
    # If in need, the scheme need to be refreshed with customized discretization info

    # Update a few dimensional parameter
    var_idxs = s[:vars]

    lbv2p = Dict()  # Bookeeping the partition idx between iterations
    for v in var_idxs
        vpcnt = length(m.discretization[v]) - 1
        chosen_p = track_new_partition_idx(m.discretization, v, m.best_bound_sol[v])
        lbv2p[v] = [i for i in 1:vpcnt if !(i in chosen_p)]
    end

    for i in 1:m.bound_sol_pool[:cnt] # First update the existing solution pool
        # First update the discretization idx with the existing
        m.bound_sol_pool[:disc][i] = Dict(j => get_active_partition_idx(m.discretization, m.bound_sol_pool[:sol][i][j],j) for j in var_idxs)
        act = true # Then check if the update pool solution active partition idex is within the deactivated region
        for v in var_idxs
            (m.bound_sol_pool[:disc][i][v] in lbv2p[v]) || (act = false)
            act || (m.bound_sol_pool[:stat][i] = :Dead)  # No re-activation
            act || break
        end
    end

    for i in 1:s[:cnt] # Now perform the merge
        act = true # Then check if the update pool solution active partition index is within the deactivated region
        for v in var_idxs
            (s[:disc][i][v] in lbv2p[v]) || (act = false)
            act || (s[:stat][i] = :Dead)
            act || break
        end
        # Reject solutions that is around best bound to avoid traps
        if isapprox(s[:obj][i], m.best_bound;atol=m.tol)
            s[:stat][i] = :Dead
        end
        push!(m.bound_sol_pool[:sol], s[:sol][i])
        push!(m.bound_sol_pool[:obj], s[:obj][i])
        push!(m.bound_sol_pool[:disc], s[:disc][i])
        push!(m.bound_sol_pool[:stat], s[:stat][i])
        push!(m.bound_sol_pool[:iter], s[:iter][i])
        push!(m.bound_sol_pool[:ubstart], s[:ubstart][i])
    end

    # Update dimensional parameters
    m.bound_sol_pool[:cnt] = length(m.bound_sol_pool[:sol])
    m.bound_sol_pool[:vars] = var_idxs

    # Show the summary
    m.loglevel > 99 && println("POOL size = $(length([i for i in 1:m.bound_sol_pool[:cnt] if m.bound_sol_pool[:stat][i] != :Dead])) / $(m.bound_sol_pool[:cnt]) ")
    for i in 1:m.bound_sol_pool[:cnt]
        m.loglevel > 99 && m.bound_sol_pool[:stat][i] != :Dead && println("ITER $(m.bound_sol_pool[:iter][i]) | SOL $(i) | POOL solution obj = $(m.bound_sol_pool[:obj][i])")
    end

    return
end

"""
    Need to be careful with the diverted point.
"""
function track_new_partition_idx(d::Dict, idx::Int, val::Float64, acp::Int=-1)

    acp > 0 ? acp = acp : acp = get_active_partition_idx(d,val,idx)

    pcnt = length(d[idx]) - 1
    newpidx = []                # Tracks the newly constructed partition idxes
    pcnt == 1 && return [1;]    # Still keep non-discretizing variables
    if acp == 1
        return newpidx = [1; 2]
    elseif acp == pcnt
        return newpidx = [pcnt-1; pcnt]
    else
        tlb = d[idx][acp-1]
        tub = d[idx][acp+1]
        if abs(val-tlb) == abs(val-tub)
            return [acp-1; acp; acp+1]
        elseif abs(val-tlb) > abs(val-tub)
            return [acp-1; acp]
        elseif abs(val-tlb) < abs(val-tub)
            return [acp; acp+1]
        end
    end

    return
end

"""
    Collect active partition idx
    Need to be careful with the diverted point
"""
function get_active_partition_idx(discretization::Dict, val::Float64, idx::Int; tol=1e-6)

    for j in 1:length(discretization[idx])-1
        if val > discretization[idx][j] - tol && val < discretization[idx][j+1] + tol
            return j
        end
    end

    @warn "Activate parition not found [VAR$(idx)]. Returning default partition 1."
    return 1
end

"""

    ncvar_collect_nodes(m:AlpineNonlinearModel)

A built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.

"""
function ncvar_collect_nodes(m::AlpineNonlinearModel;getoutput=false)

    # Pick variables that is bound width more than tolerance length
    if getoutput
        return [i for i in m.candidate_disc_vars]
    else
        m.disc_vars = [i for i in m.candidate_disc_vars]
        m.num_var_disc_mip = length(m.disc_vars)
    end

    return
end

function eval_objective(m::AlpineNonlinearModel; svec::Vector=[])

    isempty(svec) ? svec = m.best_bound_sol : svec = svec
    m.sense_orig == :Min ? obj = Inf : obj=-Inf

    if m.obj_structure == :affine
        obj = m.bounding_obj_mip[:rhs]
        for i in 1:m.bounding_obj_mip[:cnt]
            obj += m.bounding_obj_mip[:coefs][i]*svec[m.bounding_obj_mip[:vars][i].args[2]]
        end
    elseif m.structural_obj == :convex
        error("need implementation for local objective function evaluation for convex form")
    else
        error("Unknown structural obj type $(m.structural_obj)")
    end

    return obj
end

function initialize_solution_pool(m::AlpineNonlinearModel, cnt::Int)

    s = Dict()

    s[:cnt] = cnt

    # Column dimension changing variable
    s[:len] = m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip

    # !! Be careful with the :vars when utilizing the dynamic discretization variable selection !!
    s[:vars] = [i for i in m.candidate_disc_vars if length(m.discretization[i]) > 2]

    s[:sol] = Vector{Vector}(undef, cnt)                   # Solution value
    s[:obj] = Vector{Float64}(undef, cnt)                  # Objecitve value
    s[:disc] = Vector{Dict}(undef, cnt)                    # Discretization
    s[:stat] = [:Alive for i in 1:cnt]              # Solution status
    s[:iter] = [m.logs[:n_iter] for i in 1:cnt]     # Iteration collected
    s[:ubstart] = [false for i in 1:cnt]           # Solution used for ub multistart

    return s
end

"""
    Reconsideration required
"""
function ncvar_collect_arcs(m::AlpineNonlinearModel, nodes::Vector)

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
        elseif m.nonconvex_terms[k][:nonlinear_type] == :INTLIN
            var_idxs = copy(m.nonconvex_terms[k][:var_idxs])
            push!(arcs, sort(var_idxs))
        elseif m.nonconvex_terms[k][:nonlinear_type] == :INTPROD
            var_idxs = m.nonconvex_terms[k][:var_idxs]
            for i in 1:length(var_idxs)
                for j in 1:length(var_idxs)
                    i != j && push!(arcs, sort([var_idxs[i]; var_idxs[j]]))
                end
            end
        elseif m.nonconvex_terms[k][:nonlinear_type] in [:cos, :sin]
            @assert length(m.nonconvex_terms[k][:var_idxs]) == 1
            var_idx = m.nonconvex_terms[k][:var_idxs][1]
            push!(arcs, [var_idx; var_idx])
        elseif m.nonconvex_terms[k][:nonlinear_type] in [:BININT, :BINLIN, :BINPROD]
            continue
        else
            error("[EXCEPTION] Unexpected nonlinear term when building discvar graph.")
        end
    end

    return arcs
end

"""
    Special funtion for debugging bounding models
"""
function print_iis_gurobi(m::JuMP.Model)

    grb = interface_get_rawsolver(m)
    Gurobi.computeIIS(grb)
    numconstr = Gurobi.num_constrs(grb)
    numvar = Gurobi.num_vars(grb)

    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

    info("Irreducible Inconsistent Subsystem (IIS)")
    info("Variable bounds:")
    for i in 1:numvar
        v = Variable(m, i)
        if iislb[i] != 0 && iisub[i] != 0
            println(getlowerbound(v), " <= ", getname(v), " <= ", getupperbound(v))
        elseif iislb[i] != 0
            println(getname(v), " >= ", getlowerbound(v))
        elseif iisub[i] != 0
            println(getname(v), " <= ", getupperbound(v))
        end
    end

    info("Constraints:")
    for i in 1:numconstr
        if iisconstr[i] != 0
            println(m.linconstr[i])
        end
    end

    return
end

"""
    TODO can be improved
"""
function build_discvar_graph(m::AlpineNonlinearModel)

    # Collect the information of nonlinear terms in terms of arcs and nodes
    nodes = ncvar_collect_nodes(m, getoutput=true)
    arcs = ncvar_collect_arcs(m, nodes)

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

function min_vertex_cover(m::AlpineNonlinearModel)

    nodes, arcs = build_discvar_graph(m)

    # Set up minimum vertex cover problem
    update_mip_time_limit(m, timelimit=60.0)  # Set a timer to avoid waste of time in proving optimality
    minvertex = Model(solver=m.mip_solver)
    @variable(minvertex, x[nodes], Bin)
    @constraint(minvertex, [a in arcs], x[a[1]] + x[a[2]] >= 1)
    @objective(minvertex, Min, sum(x))

    status = solve(minvertex, suppress_warnings=true)
    xVal = getvalue(x)

    # Collecting required information
    m.num_var_disc_mip = Int(sum(xVal))
    m.disc_vars = [i for i in nodes if xVal[i] > m.tol && abs(m.u_var_tight[i]-m.l_var_tight[i]) >= m.tol]

    return
end

function weighted_min_vertex_cover(m::AlpineNonlinearModel, distance::Dict)

    # Collect the graph information
    nodes, arcs = build_discvar_graph(m)

    # A little bit redundency before
    disvec = [distance[i] for i in keys(distance) if i in m.candidate_disc_vars]
    disvec = abs.(disvec[disvec .> 0.0])
    isempty(disvec) ? heavy = 1.0 : heavy = 1/minimum(disvec)
    weights = Dict()
    for i in m.candidate_disc_vars
        isapprox(distance[i], 0.0; atol=1e-6) ? weights[i] = heavy : (weights[i]=(1/distance[i]))
        (m.loglevel > 100) && println("VAR$(i) WEIGHT -> $(weights[i]) ||| DISTANCE -> $(distance[i])")
    end

    # Set up minimum vertex cover problem
    update_mip_time_limit(m, timelimit=60.0)  # Set a timer to avoid waste of time in proving optimality
    minvertex = Model(solver=m.mip_solver)
    @variable(minvertex, x[nodes], Bin)
    for arc in arcs
        @constraint(minvertex, x[arc[1]] + x[arc[2]] >= 1)
    end
    @objective(minvertex, Min, sum(weights[i]*x[i] for i in nodes))

    # Solve the minimum vertex cover
    status = solve(minvertex, suppress_warnings=true)

    xVal = getvalue(x)
    m.num_var_disc_mip = Int(sum(xVal))
    m.disc_vars = [i for i in nodes if xVal[i] > 0 && abs(m.u_var_tight[i]-m.l_var_tight[i]) >= m.tol]
    m.loglevel >= 99 && println("UPDATED DISC-VAR COUNT = $(length(m.disc_vars)) : $(m.disc_vars)")

    return
end

function round_sol(m::AlpineNonlinearModel;nlp_model=nothing, nlp_sol=[])

    if nlp_model != nothing
        relaxed_sol = interface_get_solution(nlp_model)
    end

    if !isempty(nlp_sol)
        relaxed_sol = nlp_sol
    end

    if nlp_model != nothing && !isempty(nlp_sol)
        error("In function collision. Special usage")
    end

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

"""
    Evaluate a solution feasibility: Solution bust be in the feasible category and evaluated rhs must be feasible
"""
function eval_feasibility(m::AlpineNonlinearModel, sol::Vector)

    length(sol) == m.num_var_orig || error("Candidate solution length mismatch.")


    for i in 1:m.num_var_orig
        # Check solution category and bounds
        if m.var_type[i] == :Bin
            isapprox(sol[i], 1.0;atol=m.tol) || isapprox(sol[i], 0.0;atol=m.tol) || return false
        elseif m.var_type[i] == :Int
            isapprox(mod(sol[i], 1.0), 0.0;atol=m.tol) || return false
        end
        # Check solution bounds (with tight bounds)
        sol[i] <= m.l_var_tight[i] - m.tol || return false
        sol[i] >= m.u_var_tight[i] + m.tol || return false
    end

    # Check constraint violation
    eval_rhs = zeros(m.num_constr_orig)
    interface_eval_g(m.d_orig, eval_rhs, rounded_sol)
    feasible = true
    for i in 1:m.num_constr_orig
        if m.constr_type_orig[i] == :(==)
            if !isapprox(eval_rhs[i], m.l_constr_orig[i]; atol=m.tol)
                feasible = false
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) != RHS $(m.l_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        elseif m.constr_type_orig[i] == :(>=)
            if !(eval_rhs[i] >= m.l_constr_orig[i] - m.tol)
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !>= RHS $(m.l_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        elseif m.constr_type_orig[i] == :(<=)
            if !(eval_rhs[i] <= m.u_constr_orig[i] + m.tol)
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !<= RHS $(m.u_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        end
    end

    return feasible
end

function fetch_mip_solver_identifier(m::AlpineNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.mip_solver) : solverstring = override

    # Higher-level solvers: that can use sub-solvers
    if occursin("Pajarito", solverstring)
        m.mip_solver_id = "Pajarito"
        return
    elseif occursin("Pavito", solverstring)
        m.mip_solver_id = "Pavito"
        return
    end

    # Lower level solvers
    if occursin("Gurobi", solverstring)
        m.mip_solver_id = "Gurobi"
    elseif occursin("CPLEX", solverstring)
        m.mip_solver_id = "CPLEX"
    elseif occursin("Cbc", solverstring)
        m.mip_solver_id = "Cbc"
    elseif occursin("GLPK", solverstring)
        m.mip_solver_id = "GLPK"
    else
        error("Unsupported MIP solver $solverstring; use a Alpine-supported MIP solver")
    end

    return
end

function fetch_nlp_solver_identifier(m::AlpineNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.nlp_solver) : solverstring = override

    # Higher-level solver
    if occursin("Pajarito", solverstring)
        m.nlp_solver_id = "Pajarito"
        return
    elseif occursin("Pavito", solverstring)
        m.nlp_solver_id = "Pavito"
        return
    end

    # Lower-level solver
    if occursin("Ipopt", solverstring)
        m.nlp_solver_id = "Ipopt"
    elseif occursin("AmplNL", solverstring) && occursin("bonmin", solverstring)
        m.nlp_solver_id = "Bonmin"
    elseif occursin("KNITRO", solverstring)
        m.nlp_solver_id = "Knitro"
    elseif occursin("NLopt", solverstring)
        m.nlp_solver_id = "NLopt"
    else
        error("Unsupported NLP local solver $solverstring; use a Alpine-supported NLP local solver")
    end

    return
end

function fetch_minlp_solver_identifier(m::AlpineNonlinearModel;override="")

    (m.minlp_solver == empty_solver) && return

    isempty(override) ? solverstring = string(m.minlp_solver) : solverstring = override

    # Higher-level solver
    if occursin("Pajarito", solverstring)
        m.minlp_solver_id = "Pajarito"
        return
    elseif occursin("Pavito", solverstring)
        m.minlp_solver_id = "Pavito"
        return
    end

    # Lower-level Solver
    if occursin("AmplNL", solverstring) && occursin("bonmin", solverstring)
        m.minlp_solver_id = "Bonmin"
    elseif occursin("KNITRO", solverstring)
        m.minlp_solver_id = "Knitro"
    elseif occursin("NLopt", solverstring)
        m.minlp_solver_id = "NLopt"
    elseif occursin("CoinOptServices.OsilSolver(\"bonmin\"", solverstring)
        m.minlp_solver_id = "Bonmin"
    else
        error("Unsupported MINLP local solver $solverstring; use a Alpine-supported MINLP local solver")
    end

    return
end

"""
    update_mip_time_limit(m::AlpineNonlinearModel)
An utility function used to dynamically regulate MILP solver time limits to fit Alpine solver time limits.
"""
function update_mip_time_limit(m::AlpineNonlinearModel; kwargs...)

    options = Dict(kwargs)
    timelimit = 0.0
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])
    
    opts = Vector{Any}(undef, 0)
    if m.mip_solver_id != "Pavito" && m.mip_solver_id != "Pajarito"
        for i in collect(m.mip_solver.options)
            push!(opts, i)
        end
    end

    if m.mip_solver_id == "CPLEX"
        opts = update_timeleft_symbol(opts, :CPX_PARAM_TILIM, timelimit)
        m.mip_solver.options = opts
    elseif m.mip_solver_id == "Pavito"
        (timelimit < Inf) && (m.mip_solver.timeout = timelimit)
    elseif m.mip_solver_id == "Gurobi"
        opts = update_timeleft_symbol(opts, :TimeLimit, timelimit)  
        m.mip_solver.options = opts
    elseif m.mip_solver_id == "Cbc"
        opts = update_timeleft_symbol(opts, :seconds, timelimit)
        m.mip_solver.options = opts
    elseif m.mip_solver_id == "GLPK"
        opts = update_timeleft_symbol(opts, :tm_lim, timelimit)
        m.mip_solver.options = opts
    elseif m.mip_solver_id == "Pajarito"
        (timelimit < Inf) && (m.mip_solver.timeout = timelimit)
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""
    update_mip_time_limit(m::AlpineNonlinearModel)
An utility function used to dynamically regulate MILP solver time limits to fit Alpine solver time limits.
"""
function update_nlp_time_limit(m::AlpineNonlinearModel; kwargs...)

    options = Dict(kwargs)
    timelimit = 0.0
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    opts = Vector{Any}(undef, 0)
    if m.nlp_solver_id != "Pavito" && m.nlp_solver_id != "Pajarito"
        opts = collect(m.nlp_solver.options)
    end
    

    if m.nlp_solver_id == "Ipopt"
        opts = update_timeleft_symbol(opts, :max_cpu_time, timelimit)
        m.nlp_solver.options = opts
    elseif m.nlp_solver_id == "Pajarito"
        (timelimit < Inf) && (m.nlp_solver.timeout = timelimit)
    elseif m.nlp_solver_id == "AmplNL"
        opts = update_timeleft_symbol(opts, :seconds, timelimit, options_string_type=2)
        m.nlp_solver.options = opts
    elseif m.nlp_solver_id == "Knitro"
        error("You never tell me anything about knitro. Probably because they have a very short trail length.")
    elseif m.nlp_solver_id == "NLopt"
        m.nlp_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""
    update_mip_time_limit(m::AlpineNonlinearModel)
    An utility function used to dynamically regulate MILP solver time limits to fit Alpine solver time limits.
"""
function update_minlp_time_limit(m::AlpineNonlinearModel; kwargs...)

    options = Dict(kwargs)
    timelimit = 0.0
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    opts = Vector{Any}(undef, 0)
    if m.minlp_solver_id != "Pavito" && m.minlp_solver_id != "Pajarito"
        opts = collect(m.minlp_solver.options)
    end
    

    if m.minlp_solver_id == "Pajarito"
        (timelimit < Inf) && (m.minlp_solver.timeout = timelimit)
    elseif m.minlp_solver_id == "Pavito"
        (timelimit < Inf) && (m.minlp_solver.timeout = timelimit)
    elseif m.minlp_solver_id == "AmplNL"
        opts = update_timeleft_symbol(opts, :seconds, timelimit, options_string_type=2)
        m.minlp_solver.options = opts
    elseif m.minlp_solver_id == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    elseif m.minlp_solver_id == "NLopt"
        m.minlp_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""
    Follow the definition of terms to calculate the value of lifted terms
"""
function resolve_lifted_var_value(m::AlpineNonlinearModel, sol_vec::Array)

    @assert length(sol_vec) == m.num_var_orig
    sol_vec = [sol_vec; fill(NaN, m.num_var_linear_mip+m.num_var_nonlinear_mip)]

    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            lvar_idx = m.nonconvex_terms[k][:y_idx]
            sol_vec[lvar_idx] = m.nonconvex_terms[k][:evaluator](m.nonconvex_terms[k], sol_vec)
        elseif haskey(m.linear_terms, k)
            lvar_idx = m.linear_terms[k][:y_idx]
            sol_vec[lvar_idx] = m.linear_terms[k][:evaluator](m.linear_terms[k], sol_vec)
        else
            error("[RARE] Found homeless term key $(k) during bound resolution.")
        end
    end

    return sol_vec
end

function adjust_branch_priority(m::AlpineNonlinearModel)

    if m.mip_solver_id == "Gurobi"
        !m.model_mip.internalModelLoaded && return
        len = length(m.model_mip.colVal)
        prior = Cint[] # priorities
        for i=1:len
            push!(prior, i)
        end
        Gurobi.set_intattrarray!(m.model_mip.internalModel.inner, "BranchPriority", 1, len, prior)
    elseif m.mip_solver_id == "CPLEX"
        !m.model_mip.internalModelLoaded && return
        n = length(m.model_mip.colVal)
        idxlist = Cint[1:n;] # variable indices
        prior = Cint[] # priorities
        for i=1:n
            push!(prior, i)
        end
        CPLEX.set_branching_priority(MathProgBase.getrawsolver(internalmodel(m.model_mip)), idxlist, prior)
    else
        return
    end

end
