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

    if m.best_obj in [Inf, -Inf]
        m.best_rel_gap = Inf
        return
    else
        p = round(abs(log(10,m.relgap)))
        n = round(abs(m.best_obj-m.best_bound), Int(p))
        dn = round(abs(1e-12+abs(m.best_obj)), Int(p))
        if (n == 0.0) && (dn == 0.0)
            m.best_rel_gap = 0.0
            return
        end
        m.best_rel_gap = abs(m.best_obj - m.best_bound)/(m.tol+abs(m.best_obj))
    end

    # absoluate or any other bound calculation shows here..
    return
end

"""
    discretization_to_bounds(d::Dict, l::Int)

Same as [`update_var_bounds`](@ref)
"""
discretization_to_bounds(d::Dict, l::Int) = update_var_bounds(d, len=l)


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
function insert_timeleft_symbol(options, val::Float64, keywords::Symbol, timeout; options_string_type=1)

    for i in 1:length(options)
        if options_string_type == 1
            if keywords in collect(options[i])
                deleteat!(options, i)
            end
        elseif options_string_type == 2
            if keywords == split(options[i],"=")[1]
                deleteat!(options, i)
            end
        end
    end

    if options_string_type == 1
        (val != Inf) && push!(options, (keywords, val))
    elseif options_string_type == 2
        (val != Inf) && push!(options, "$(keywords)=$(val)")
    end

    return
end

"""
    fetch_boundstop_symbol(m::PODNonlinearModel)

An utility function used to recongize different sub-solvers and return the bound stop option key words
"""
function update_boundstop_options(m::PODNonlinearModel)

    if m.mip_solver_id == "Gurobi"
        # Calculation of the bound
        if m.sense_orig == :Min
            stopbound = (1-m.relgap+m.tol) * m.best_obj
        elseif m.sense_orig == :Max
            stopbound = (1+m.relgap-m.tol) * m.best_obj
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
    check_solution_history(m::PODNonlinearModel, ind::Int)

Check if the solution is alwasy the same within the last disc_consecutive_forbid iterations. Return true if suolution in invariant.
"""
function check_solution_history(m::PODNonlinearModel, ind::Int)

    (m.logs[:n_iter] < m.disc_consecutive_forbid) && return false

    sol_val = m.bound_sol_history[mod(m.logs[:n_iter]-1, m.disc_consecutive_forbid)+1][ind]
    for i in 1:(m.disc_consecutive_forbid-1)
        search_pos = mod(m.logs[:n_iter]-1-i, m.disc_consecutive_forbid)+1
        !isapprox(sol_val, m.bound_sol_history[search_pos][ind]; atol=m.disc_rel_width_tol) && return false
    end
    return true
end

"""

    fix_domains(m::PODNonlinearModel)

This function is used to fix variables to certain domains during the local solve process in the [`global_solve`](@ref).
More specifically, it is used in [`local_solve`](@ref) to fix binary and integer variables to lower bound solutions
and discretizing varibles to the active domain according to lower bound solution.
"""
function fix_domains(m::PODNonlinearModel;discrete_sol=nothing)

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
                    l_var[i] = m.discretization[i][j]
                    u_var[i] = m.discretization[i][j+1]
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
    is_fully_convexified(m::PODNonlinearModel)
"""
function is_fully_convexified(m::PODNonlinearModel)

    # Other more advanced convexification check goes here
    for term in keys(m.nonconvex_terms)
        if !m.nonconvex_terms[term][:convexified]
            warn("Detected terms that is not convexified $(term[:lifted_constr_ref]), bounding model solver may report a error due to this")
            return
        else
            m.nonconvex_terms[term][:convexified] = false    # Reset status for next iteration
        end
    end

    return
end

"""

    ncvar_collect_nodes(m:PODNonlinearModel)

A built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.

"""
function ncvar_collect_nodes(m::PODNonlinearModel;getoutput=false)

    # Pick variables that is bound width more than tolerance length
    if getoutput
        return [i for i in m.candidate_disc_vars]
    else
        m.disc_vars = [i for i in m.candidate_disc_vars]
        m.num_var_disc_mip = length(m.disc_vars)
    end

    return
end

"""
    Reconsideration required
"""
function ncvar_collect_arcs(m::PODNonlinearModel, nodes::Vector)

    arcs = Set()

    for k in keys(m.nonconvex_terms)
        if m.nonconvex_terms[k][:nonlinear_type] == :BILINEAR
            arc = [i.args[2] for i in k]
            length(arc) == 2 && push!(arcs, sort(arc))
        elseif m.nonconvex_terms[k][:nonlinear_type] == :MONOMIAL
            @assert isa(m.nonconvex_terms[k][:var_idxs][1], Int)
            varidx = m.nonconvex_terms[k][:var_idxs][1]
            push!(arcs, [varidx, varidx;])
        elseif m.nonconvex_terms[k][:nonlinear_type] == :MULTILINEAR
            varidxs = m.nonconvex_terms[k][:var_idxs]
            for i in 1:length(varidxs)
                for j in 1:length(varidxs)
                    if i != j
                        push!(arcs, sort([varidxs[i], varidxs[j];]))
                    end
                end
            end
        elseif m.nonconvex_terms[k][:nonlinear_type] == :INTLIN
            var_idxs = copy(m.nonconvex_terms[k][:var_idxs])
            push!(arcs, sort(var_idxs))
        elseif m.nonconvex_terms[k][:nonlinear_type] == :INTPROD
            var_idxs = m.nonconvex_terms[k][:var_idxs]
            for i in 1:length(var_idxs)
                for j in 1:length(var_idxs)
                    i != j && push!(arcs, sort([var_idxs[i], var_idxs[j];]))
                end
            end
        elseif m.nonconvex_terms[k][:nonlinear_type] in [:cos, :sin]
            @assert length(m.nonconvex_terms[k][:var_idxs]) == 1
            var_idx = m.nonconvex_terms[k][:var_idxs][1]
            push!(arcs, [var_idx, var_idx;])
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

    grb = MathProgBase.getrawsolver(internalmodel(m))
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
function build_discvar_graph(m::PODNonlinearModel)

    # Collect the information of nonlinear terms in terms of arcs and nodes
    nodes = ncvar_collect_nodes(m, getoutput=true)
    arcs = ncvar_collect_arcs(m, nodes)

    # Collect integer variables
    for i in 1:m.num_var_orig
        if !(i in nodes) && m.var_type[i] == :Int
            push!(nodes, i)
            push!(arcs, [i,i;])
        end
    end

    nodes = collect(nodes)
    arcs = collect(arcs)

    return nodes, arcs
end

function min_vertex_cover(m::PODNonlinearModel)

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

function weighted_min_vertex_cover(m::PODNonlinearModel, distance::Dict)

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

function round_sol(m::PODNonlinearModel, nlp_model)

    relaxed_sol = MathProgBase.getsolution(nlp_model)
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
function eval_feasibility(m::PODNonlinearModel, sol::Vector)

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
    MathProgBase.eval_g(m.d_orig, eval_rhs, rounded_sol)
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

function fetch_mip_solver_identifier(m::PODNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.mip_solver) : solverstring=override

    # Higher-level solvers: that can use sub-solvers
    if contains(solverstring,"Pajarito")
        m.mip_solver_id = "Pajarito"
        return
    end

    # Lower level solvers
    if contains(solverstring,"Gurobi")
        m.mip_solver_id = "Gurobi"
    elseif contains(solverstring,"CPLEX")
        m.mip_solver_id = "CPLEX"
    elseif contains(solverstring,"Cbc")
        m.mip_solver_id = "Cbc"
    elseif contains(solverstring,"GLPK")
        m.mip_solver_id = "GLPK"
    else
        error("Unsupported mip solver name. Using blank")
    end

    return
end

function fetch_nlp_solver_identifier(m::PODNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.nlp_solver) : solverstring=override

    # Higher-level solver
    if contains(solverstring, "Pajarito")
        m.nlp_solver_id = "Pajarito"
        return
    end

    # Lower-level solver
    if contains(solverstring, "Ipopt")
        m.nlp_solver_id = "Ipopt"
    elseif contains(solverstring, "AmplNL") && contains(solverstring, "bonmin")
        m.nlp_solver_id = "Bonmin"
    elseif contains(solverstring, "KNITRO")
        m.nlp_solver_id = "Knitro"
    elseif contains(solverstring, "NLopt")
        m.nlp_solver_id = "NLopt"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

function fetch_minlp_solver_identifier(m::PODNonlinearModel;override="")

    (m.minlp_solver == UnsetSolver()) && return

    isempty(override) ? solverstring = string(m.minlp_solver) : solverstring=override

    # Higher-level solver
    if contains(solverstring, "Pajarito")
        m.minlp_solver_id = "Pajarito"
        return
    end

    # Lower-level Solver
    if contains(solverstring, "AmplNL") && contains(solverstring, "bonmin")
        m.minlp_solver_id = "Bonmin"
    elseif contains(solverstring, "KNITRO")
        m.minlp_solver_id = "Knitro"
    elseif contains(solverstring, "NLopt")
        m.minlp_solver_id = "NLopt"
    elseif contains(solverstring, "CoinOptServices.OsilSolver(\"bonmin\"")
        m.minlp_solver_id = "Bonmin"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_mip_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.mip_solver_id == "CPLEX"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.mip_solver_id == "Gurobi"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:TimeLimit,m.timeout)
    elseif m.mip_solver_id == "Cbc"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:seconds,m.timeout)
    elseif m.mip_solver_id == "GLPK"
        insert_timeleft_symbol(m.mip_solver.opts, timelimit,:tm_lim,m.timeout)
    elseif m.mip_solver_id == "Pajarito"
        (timelimit < Inf) && (m.mip_solver.timeout = timelimit)
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_nlp_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.nlp_solver_id == "Ipopt"
        insert_timeleft_symbol(m.nlp_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.nlp_solver_id == "Pajarito"
        (timelimit < Inf) && (m.nlp_solver.timeout = timelimit)
    elseif m.nlp_solver_id == "AmplNL"
        insert_timeleft_symbol(m.nlp_solver.options,timelimit,:seconds,m.timeout, options_string_type=2)
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

    update_mip_time_limit(m::PODNonlinearModel)

    An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_minlp_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.minlp_solver_id == "Pajarito"
        (timelimit < Inf) && (m.minlp_solver.timeout = timelimit)
    elseif m.minlp_solver_id == "AmplNL"
        insert_timeleft_symbol(m.minlp_solver.options,timelimit,:seconds,m.timeout,options_string_type=2)
    elseif m.minlp_solver_id == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    elseif m.minlp_solver_id == "NLopt"
        m.minlp_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end
