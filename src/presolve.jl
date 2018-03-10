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

    if m.presolve_bt == false # User choose not to do bound tightening
        return
    end

    if m.presolve_bt_algo == 1
        minmax_bound_tightening(m, use_bound=use_bound)
    elseif m.presolve_bt_algo == 2
        minmax_bound_tightening(m, use_bound=use_bound, use_tmc=true)
    elseif isa(m.presolve_bt_algo, Function)
        eval(m.presolve_bt_algo)(m)
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
    status_pass = [:Optimal]
    status_reroute = [:UserLimits]

    options = Dict(kwargs)

    # Regulating Special Input Conditions: default use best feasible solution objective value
    (use_bound == true) ? bound = m.best_obj : bound = Inf
    l_var_orig = copy(m.l_var_tight)
    u_var_orig = copy(m.u_var_tight)

    discretization = to_discretization(m, m.l_var_tight, m.u_var_tight)
    if use_bound == false && haskey(options, :use_tmc)
        (m.loglevel > 0) && warn("[BOUND TIGHTENING ALGO] TMC chosen by the user, but local solve infeasible; defaulting to doing bound-tightening without TMC.")
    end
    if use_bound == true && haskey(options, :use_tmc)
        discretization = add_adaptive_partition(m, use_solution=m.best_sol, use_disc=discretization)
    end
    discretization = resolve_var_bounds(m, discretization) # recomputation of bounds for lifted_variables

    (m.loglevel > 0) && println("starting the bound-tightening algorithm ...")

    # start of the solve
    keeptightening = true
    while keeptightening && (m.logs[:time_left] > m.tol) && (m.logs[:bt_iter] < m.presolve_maxiter) # Stopping criteria

        keeptightening = false
        m.logs[:bt_iter] += 1
        m.loglevel > 99 && println("[DEBUG] Iteration - $(m.logs[:bt_iter])")
        temp_bounds = Dict()

        # Perform Bound Contraction
        for var_idx in m.candidate_disc_vars # why only all nonlinear vars?
            temp_bounds[var_idx] = [discretization[var_idx][1], discretization[var_idx][end]]
            if abs(discretization[var_idx][1] - discretization[var_idx][end]) > m.presolve_bt_width_tol
                create_bound_tightening_model(m, discretization, bound)
                for sense in both_senses
                    @objective(m.model_mip, sense, Variable(m.model_mip, var_idx))
                    status = solve_bound_tightening_model(m)
                    if status in status_pass
                        temp_bounds[var_idx][tell_side[sense]] = eval(tell_round[sense])(getobjectivevalue(m.model_mip)/m.presolve_bt_output_tol)*m.presolve_bt_output_tol  # Objective truncation for numerical issues
                    elseif status in status_reroute
                        temp_bounds[var_idx][tell_side[sense]] = eval(tell_round[sense])(getobjbound(m.model_mip)/m.presolve_bt_output_tol)*m.presolve_bt_output_tol
                    else
                        print("!")
                        temp_bounds[var_idx][tell_side[sense]] = temp_bounds[var_idx][tell_side[sense]]
                    end
                    m.loglevel > 199 && println("[DEBUG] contracting VAR $(var_idx) $(sense) problem, results in $(temp_bounds[var_idx][tell_side[sense]])")
                end
            end
        end

        # Updates the discretization structure
        for var_idx in m.candidate_disc_vars
            if temp_bounds[var_idx][1] >= discretization[var_idx][1] + m.presolve_bt_width_tol
                (m.loglevel > 99) && print("+")
                keeptightening = true # Continue to perform the next iteration
                discretization[var_idx][1] = temp_bounds[var_idx][1]
            end
            if discretization[var_idx][end] <= temp_bounds[var_idx][end] - m.presolve_bt_width_tol
                (m.loglevel > 99) && print("+")
                keeptightening = true
                discretization[var_idx][end] = temp_bounds[var_idx][end]
            end
        end
        (m.loglevel > 0) && print("\n")
        discretization = resolve_var_bounds(m, discretization)
        haskey(options, :use_tmc) ? discretization = add_adaptive_partition(m, use_solution=m.best_sol, use_disc=flatten_discretization(discretization)) : discretization = discretization
    end

    (m.loglevel > 0) && println("\nfinished bound tightening in $(m.logs[:bt_iter]) iterations, applying tighten bounds")

    m.l_var_tight, m.u_var_tight = update_var_bounds(discretization)
    m.discretization = add_adaptive_partition(m, use_solution=m.best_sol)

    (m.loglevel > 0)  && [println("[DEBUG] VAR $(i) BOUND contracted $(round(1-abs(m.l_var_tight[i] - m.u_var_tight[i])/abs(l_var_orig[i] - u_var_orig[i]),2)*100)% |$(round(l_var_orig[i],4)) --> | $(round(m.l_var_tight[i],4)) - * - $(round(m.u_var_tight[i],4)) | <-- $(round(u_var_orig[i],4)) |") for i in 1:m.num_var_orig]
    (m.loglevel > 0) && print("\n")
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
    amp_post_vars(m, use_disc=discretization)
    amp_post_lifted_constraints(m)
    amp_post_convexification(m, use_disc=discretization)  # Convexify problem

    if bound != Inf
        post_obj_bounds(m, bound)
    end

    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

    return
end

"""

    solve_bound_tightening_model(m::PODNonlinearModel)

A function that solves the min and max bound-tightening model.

"""
function solve_bound_tightening_model(m::PODNonlinearModel; kwargs...)

    # ========= MILP Solve ========= #
    if m.presolve_bt_mip_timeout < Inf
        update_mip_time_limit(m, timelimit = max(0.0, min(m.presolve_bt_mip_timeout, m.timeout - m.logs[:total_time])))
    else
        update_mip_time_limit(m, timelimit = max(0.0, m.timeout - m.logs[:total_time]))
    end

    start_solve = time()
    status = solve(m.model_mip, suppress_warnings=true, relaxation=m.presolve_bt_relax) #TODO Double check here
    cputime_solve = time() - start_solve
    m.logs[:total_time] += cputime_solve
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])
    # ========= MILP Solve ========= #

    return status
end

"""
    TODO: docstring
"""
function post_obj_bounds(m::PODNonlinearModel, bound::Float64; kwargs...)
    if m.sense_orig == :Max
        @constraint(m.model_mip,
            sum(m.bounding_obj_mip[:coefs][j]*Variable(m.model_mip, m.bounding_obj_mip[:vars][j].args[2]) for j in 1:m.bounding_obj_mip[:cnt]) >= bound)
    elseif m.sense_orig == :Min
        @constraint(m.model_mip,
            sum(m.bounding_obj_mip[:coefs][j]*Variable(m.model_mip, m.bounding_obj_mip[:vars][j].args[2]) for j in 1:m.bounding_obj_mip[:cnt]) <= bound)
    end

    return
end
