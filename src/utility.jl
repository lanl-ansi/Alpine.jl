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

    @assert length(lbs) == length(ubs)
    var_discretization = Dict()
    for var in 1:m.num_var_orig
        lb = lbs[var]
        ub = ubs[var]
        var_discretization[var] = [lb, ub]
    end

    if length(lbs) == (m.num_var_orig+m.num_var_lifted_mip)
        for var in (1+m.num_var_orig):(m.num_var_orig+m.num_var_lifted_mip)
            lb = lbs[var]
            ub = ubs[var]
            var_discretization[var] = [lb, ub]
        end
    else
        for var in (1+m.num_var_orig):(m.num_var_orig+m.num_var_lifted_mip)
            lb = -Inf
            ub = Inf
            var_discretization[var] = [lb, ub]
        end
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

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_mip_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.mip_solver_identifier == "CPLEX"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.mip_solver_identifier == "Gurobi"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:TimeLimit,m.timeout)
    elseif m.mip_solver_identifier == "Cbc"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:seconds,m.timeout)
    elseif m.mip_solver_identifier == "GLPK"
        insert_timeleft_symbol(m.mip_solver.opt)
    else found == nothing
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

    if m.nlp_local_solver_identifier == "Ipopt"
        insert_timeleft_symbol(m.nlp_local_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.nlp_local_solver_identifier == "Pajarito"
        (timeout < Inf) && (m.nlp_local_solver.timeout = timelimit)
    elseif m.nlp_local_solver_identifier == "AmplNL"
        insert_timeleft_symbol(m.nlp_local_solver.options,timelimit,:seconds,m.timeout, options_string_type=2)
    elseif m.nlp_local_solver_identifier == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    else found == nothing
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

    if m.minlp_local_solver_identifier == "Pajarito"
        (timeout < Inf) && (m.minlp_local_solver.timeout = timelimit)
    elseif m.minlp_local_solver_identifier == "AmplNL"
        insert_timeleft_symbol(m.minlp_local_solver.options,timelimit,:seconds,m.timeout,options_string_type=2)
    elseif m.minlp_local_solver_identifier == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    else found == nothing
        error("Needs support for this MIP solver")
    end

    return
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
        (timeout != Inf) && push!(options, (keywords, val))
    elseif options_string_type == 2
        (timeout != Inf) && push!(options, "$(keywords)=$(val)")
    end
    return
end

"""
    fetch_boundstop_symbol(m::PODNonlinearModel)

An utility function used to recongize different sub-solvers and return the bound stop option key words
"""
function update_boundstop_options(m::PODNonlinearModel)


    # # Calculation of the bound
    # if m.sense_orig == :Min
    #     stopbound = (1-m.rel_gap+m.tol) * m.best_obj
    # elseif m.sense_orig == :Max
    #     stopbound = (1+m.rel_gap-m.tol) * m.best_obj
    # end
    #
    # for i in 1:length(m.mip_solver.options)
    #     if m.mip_solver.options[i][1] == :BestBdStop
    #         deleteat!(m.mip_solver.options, i)
    #         if string(m.mip_solver)[1:6] == "Gurobi"
    #             push!(m.mip_solver.options, (:BestBdStop, stopbound))
    #         else
    #             return
    #         end
    #     end
    # end
    #
    # if string(m.mip_solver)[1:6] == "Gurobi"
    #     push!(m.mip_solver.options, (:BestBdStop, stopbound))
    # else
    #     return
    # end

    return
end


"""
    check_solution_history(m::PODNonlinearModel, ind::Int)

Check if the solution is alwasy the same within the last discretization_consecutive_forbid iterations. Return true if suolution in invariant.
"""
function check_solution_history(m::PODNonlinearModel, ind::Int)

    (m.logs[:n_iter] < m.discretization_consecutive_forbid) && return false

    sol_val = m.bound_sol_history[mod(m.logs[:n_iter]-1, m.discretization_consecutive_forbid)+1][ind]
    for i in 1:(m.discretization_consecutive_forbid-1)
        search_pos = mod(m.logs[:n_iter]-1-i, m.discretization_consecutive_forbid)+1
        !isapprox(sol_val, m.bound_sol_history[search_pos][ind]; atol=m.discretization_rel_width_tol) && return false
    end
    return true
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
            error("Detected terms that is not convexified $(term[:lifted_constr_ref]), bounding model solver may report a error due to this")
            return
        else
            m.nonlinear_info[term][:convexified] = false    # Reset status for next iteration
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

    return
end

"""

    max_cover(m:PODNonlinearModel)

A built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.

"""
function max_cover(m::PODNonlinearModel; kwargs...)

    nodes = Set()
    for pair in keys(m.nonlinear_info)
        # Assumption Max cover is always safe
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

function fetch_mip_solver_identifier(m::PODNonlinearModel)

    @show string(m.nlp_local_solver)
    if string(m.mip_solver)[1:6] == "Gurobi"
        m.mip_solver_identifier = "Gurobi"
    elseif string(m.mip_solver)[1:7] == "CPLEX"
        m.mip_solver_identifier = "Cplex"
    elseif string(m.mip_solver)[1:3] == "Cbc"
        m.mip_solver_identifier = "Cbc"
    elseif string(m.mip_solver)[1:4] == "GLPK"
        m.mip_solver_identifier = "GLPK"
    else
        error("Unsupported mip solver name. Using blank")
    end

    return
end

function fetch_nlp_solver_identifier(m::PODNonlinearModel)

    @show string(m.nlp_local_solver)
    if string(m.nlp_local_solver)[1:5] == "Ipopt"
        m.nlp_local_solver_identifier = "Ipopt"
    elseif string(m.nlp_local_solver)[1:6] == "AmplNL"
        m.nlp_local_solver_identifier = "Bonmin"
    elseif string(m.nlp_local_solver)[1:6] == "Knitro"
        m.nlp_local_solver_identifier = "Knitro"
    elseif string(m.nlp_local_solver)[1:8] == "Pajarito"
        m.nlp_local_solver_identifier = "Pajarito"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

function fetch_minlp_solver_identifier(m::PODNonlinearModel)

    (m.minlp_local_solver == UnsetSolver()) && return
    @show string(m.minlp_local_solver)
    if string(m.minlp_local_solver)[1:6] == "AmplNL"
        m.minlp_local_solver_identifier = "Bonmin"
    elseif string(m.minlp_local_solver)[1:6] == "Knitro"
        m.minlp_local_solver_identifier = "Knitro"
    elseif string(m.minlp_local_solver)[1:8] == "Pajarito"
        m.minlp_local_solver_identifier = "Pajarito"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end


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
