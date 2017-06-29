"""

    bound_tightening(m::PODNonlinearModel)

Entry point for the bound-tightening algorithm. The aim of the bound-tightening algorithm
is to tighten the variable bounds, if possible.

Currently, two bounding tightening method is implemented [`minmax_bound_tightening`](@ref).

    * Bound-tightening with basic McCormick
    * Bound-tightening with McCormick partitions: (3 partitions around the local feasible solution)
    If no local feasible solution is obtained, the algorithm defaults to bound-tightening with basic McCormick

"""
function bound_tightening(m::PODNonlinearModel; use_bound = true, kwargs...)

    if m.presolve_perform_bound_tightening == false # User choose not to do bound tightening
        return
    end

    if m.presolve_bound_tightening_algo == 1
        minmax_bound_tightening(m, use_bound=use_bound)
    elseif m.presolve_bound_tightening_algo == 2
        minmax_bound_tightening(m, use_bound=use_bound, use_tmc=true)
    elseif isa(m.presolve_bound_tightening_algo, Function)
        eval(m.presolve_bound_tightening_algo)(m)
    else
        error("Unrecognized bound-tightening algorithm")
    end

    return
end

"""

    minmax_bound_tightening(m:PODNonlinearModel; use_bound::Bool=true, use_tmc::Bool)

This function implements the bound-tightening algorithm to tighten the variable bounds.
It utilizes either the basic McCormick relaxation or the Tightened McCormick relaxation (TMC)
to tighten the bounds. The TMC has additional binary variables for partitioning.

The algorithm as two main parameters. The first is the `use_tmc`, which when set to `true`
invokes the algorithm on the TMC relaxation. The second parameter `use_bound` takes in the
objective value of the local solve solution stored in `best_sol`. The `use_bound` option is set
to `true` when the local solve is successful is obtaining a feasible solution and the bound-tightening
is to be performed using the objective value of the feasible solution, else this parameter
is set to `false`

Several other parameters are available for the presolve algorithm tuning.
For more details, see [Parameters](@ref).

"""
function minmax_bound_tightening(m::PODNonlinearModel; use_bound = true, kwargs...)

    # Some functinal constants
    both_senses = [:Min, :Max]             # Senses during bound tightening procedures
    tell_side = Dict(:Min=>1, :Max=>2)     # Positional information
    tell_round = Dict(:Min=>floor, :Max=>ceil)

    options = Dict(kwargs)

    # Regulating Special Input Conditions: default use best feasible solution objective value
    (use_bound == true) ? bound = m.best_obj : bound = Inf
    discretization = to_discretization(m, m.l_var_tight, m.u_var_tight)
    if use_bound == false && haskey(options, :use_tmc)
        (m.log_level > 0) && warn("[BOUND TIGHTENING ALGO] TMC chosen by the user, but local solve infeasible; defaulting to doing bound-tightening without TMC.")
    end
    if use_bound == true && haskey(options, :use_tmc)
        discretization = add_discretization(m, use_solution=m.best_sol, use_discretization=discretization)
    end
    discretization = resolve_lifted_var_bounds(m.nonlinear_info, discretization) # recomputation of bounds for lifted_variables


    (m.log_level > 0) && println("starting the bound-tightening algorithm ...")
    (m.log_level > 99) && [println("[DEBUG] VAR $(var_idx) Original Bound [$(round(m.l_var_tight[var_idx],4)) < - > $(round(m.u_var_tight[var_idx],4))]") for var_idx in m.all_nonlinear_vars]

    # start of the solve
    keeptightening = true
    while keeptightening && m.logs[:time_left] > m.tolerance && m.logs[:bt_iter] < m.presolve_maxiter # Stopping criteria

        keeptightening = false
        m.logs[:bt_iter] += 1
        (m.log_level > 99) && println("[DEBUG] Iteration - $(m.logs[:bt_iter])")
        temp_bounds = Dict()

        # Perform Bound Contraction
        for var_idx in m.all_nonlinear_vars # why only all nonlinear vars?
            temp_bounds[var_idx] = [discretization[var_idx][1], discretization[var_idx][end]]
            if abs(discretization[var_idx][1] - discretization[var_idx][end]) > m.presolve_bt_width_tolerance
                create_bound_tightening_model(m, discretization, bound)
                for sense in both_senses
                    @objective(m.model_mip, sense, Variable(m.model_mip, var_idx))
                    solve_bound_tightening_model(m)
                    temp_bounds[var_idx][tell_side[sense]] = eval(tell_round[sense])(getobjectivevalue(m.model_mip)/m.presolve_bt_output_tolerance)*m.presolve_bt_output_tolerance  # Objective truncation for numerical issues
                    # TODO: discuss feasibility tolerances and where to put them and apt default
                    m.log_level > 99 && println("[DEBUG] contracting VAR $(var_idx) $(sense) problem, results in $(temp_bounds[var_idx][tell_side[sense]])")
                end
            end
        end

        # Updates the discretization structure
        for var_idx in m.all_nonlinear_vars
            if abs((temp_bounds[var_idx][1] - discretization[var_idx][1])/discretization[var_idx][1]) > m.presolve_bt_width_tolerance
                keeptightening = true # Continue to perform the next iteration
                discretization[var_idx][1] = temp_bounds[var_idx][1]
            end
            if abs((discretization[var_idx][end]-temp_bounds[var_idx][end])/discretization[var_idx][end]) > m.presolve_bt_width_tolerance
                (m.log_level > 0) && print("+")
                keeptightening = true
                discretization[var_idx][end] = temp_bounds[var_idx][end]
            end
        end

        discretization = resolve_lifted_var_bounds(m.nonlinear_info, discretization)
        haskey(options, :use_tmc) ? discretization = add_discretization(m, use_solution=m.best_sol, use_discretization=flatten_discretization(discretization)) : discretization = discretization
    end

    (m.log_level > 0) && println("\nfinished bound tightening in $(m.logs[:bt_iter]) iterations, applying tighten bounds")

    # Update the bounds with the tightened ones
    m.l_var_tight, m.u_var_tight = update_var_bounds(discretization)

    (m.log_level > 99)  && [println("[DEBUG] VAR $(i) BOUND contracted |$(round(m.l_var_orig[i],4)) --> | $(round(m.l_var_tight[i],4)) - * - $(round(m.u_var_tight[i],4)) | <-- $(round(m.u_var_orig[i],4)) |") for i in m.all_nonlinear_vars]
    (m.log_level > 0) && print("\n")
    return
end

"""
    create_bound_tightening_model(m::PODNonlinearModel, discretization::Dict, bound::Float64)

This function takes in the initial discretization information and builds a bound-tightening model.
It is an algorithm specific function called by [`minmax_bound_tightening`](@ref)

 """
function create_bound_tightening_model(m::PODNonlinearModel, discretization, bound; kwargs...)

    options = Dict(kwargs)

    start_build = time()

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP model
    post_amp_vars(m, use_discretization=discretization)
    post_amp_lifted_constraints(m)
    post_amp_mccormick(m, use_discretization=discretization)
    # any additional built-in convexification method : convexhull for example
    for i in 1:length(m.expression_convexification)             # Additional user-defined convexificaition method
        eval(m.expression_convexification[i])(m)
    end

    if bound != Inf
        post_obj_bounds(m, bound)
    end

    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build * m.presolve_track_time
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time] * m.presolve_track_time)

    return
end

"""

    solve_bound_tightening_model(m::PODNonlinearModels)

A function that solves the min and max bound-tightening model.

"""
function solve_bound_tightening_model(m::PODNonlinearModel; kwargs...)

    # ========= MILP Solve ========= #
    if m.presolve_mip_timelimit < Inf
        update_mip_time_limit(m, timelimit = max(0.0, min(m.presolve_mip_timelimit, m.timeout - m.logs[:total_time])))
    else
        update_mip_time_limit(m, timelimit = max(0.0, m.timeout - m.logs[:total_time]))
    end

    start_solve = time()
    status = solve(m.model_mip, suppress_warnings=true, relaxation=m.presolve_mip_relaxation) #TODO Double check here
    cputime_solve = time() - start_solve
    m.logs[:total_time] += cputime_solve * m.presolve_track_time
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time] * m.presolve_track_time)
    #TODO handle the infeasible cases when time limit is applied
    # ========= MILP Solve ========= #

    return
end

"""
    resolve_lifted_var_bounds(nonlinear_info::Dict, discretization::Dict)

For discretization to be performed, we do not allow for a variable being discretized to have infinite bounds.
The lifted variables will have infinite bounds and the function infers bounds on these variables. This process
can help speed up the subsequent solve in subsequent iterations.
"""
function resolve_lifted_var_bounds(nonlinear_info::Dict, discretization::Dict; kwargs...)

    options = Dict(kwargs)

    for bi in keys(nonlinear_info)
        idx_a = bi[1].args[2]
        idx_b = bi[2].args[2]
        idx_ab = nonlinear_info[bi][:lifted_var_ref].args[2]
        bound = [discretization[idx_a][1], discretization[idx_a][end]] * [discretization[idx_b][1], discretization[idx_b][end]]'
        discretization[idx_ab] = [-Inf, Inf]
        discretization[idx_ab][1] = minimum(bound)
        discretization[idx_ab][2] = maximum(bound)
    end

    return discretization
end


"""
    resolve_closed_var_bounds(m::PODNonlinearModel)

This function seeks variable with tight bounds (by presolve_bt_width_tolerance) by checking .l_var_tight and .u_var_tight.
If a variable is found to be within a sufficiently small interval then no discretization will be performed on this variable
and the .discretization will be cleared with the tight bounds for basic McCormick operation if necessary.

"""
function resolve_closed_var_bounds(m::PODNonlinearModel; kwargs...)

    for var in m.all_nonlinear_vars
        if abs(m.l_var_tight[var] - m.u_var_tight[var]) < m.presolve_bt_width_tolerance   # Closed Bound Criteria
            deleteat!(m.var_discretization_mip, findfirst(m.var_discretization_mip, var)) # Clean nonlinear_info by deleting the info
            m.discretization[var] = [m.l_var_tight[var], m.u_var_tight[var]]              # Clean up the discretization for basic McCormick if necessary
        end
    end

    return
end
