"""
    High-level Function
"""
function MathProgBase.optimize!(m::PODNonlinearModel)
    any(isnan, m.best_sol) && (m.best_sol = zeros(length(m.best_sol)))
    presolve(m)
    global_solve(m)
    m.loglevel > 0 && logging_row_entry(m, finsih_entry=true)
    summary_status(m)
end

"""
    global_solve(m::PODNonlinearModel)

Perform the global algorithm that is based on the adaptive conexification scheme.
This iterative algorithm loops over [`bounding_solve`](@ref) and [`local_solve`](@ref) for converging lower bound (relaxed problem) and upper bound (feasible problem).
Each [`bounding_solve`](@ref) provides a lower bound solution that is used as a partioning point for next iteration (this feature can be modified given different `add_adaptive_partition`).
Each [`local_solve`](@ref) provides a local serach of incumbent feasible solution. The algrithm terminates given time limits, optimality condition, or iteration limits.

The algorithm is can be reformed when `add_adaptive_partition` is replaced with user-defined functional input.
For example, this algorithm can easily be reformed as a uniform-partitioning algorithm in other literature.

"""
function global_solve(m::PODNonlinearModel)

    m.loglevel > 0 && logging_head(m)
    m.presolve_track_time || reset_timer(m)

    while !check_exit(m)
        m.logs[:n_iter] += 1
        create_bounding_mip(m)                  # Build the relaxation model
        bounding_solve(m)                       # Solve the relaxation model
        update_opt_gap(m)                       # Update optimality gap
        check_exit(m) && break                  # Feasibility check
        m.loglevel > 0 && logging_row_entry(m)  # Logging
        local_solve(m)                          # Solve local model for feasible solution
        update_opt_gap(m)                       # Update optimality gap
        check_exit(m) && break                  # Detect optimality termination
        algorithm_automation(m)                 # Automated adjustments
        add_partition(m)                        # Add extra discretizations
    end

    return
end


"""
    presolve(m::PODNonlinearModel)
"""
function presolve(m::PODNonlinearModel)

    start_presolve = time()
    (m.loglevel > 0) && println("\nPOD algorithm presolver started.")
    (m.loglevel > 0) && println("performing local solve to obtain a feasible solution.")
    local_solve(m, presolve = true)

    # Possible solver status, return error when see different
    status_pass = [:Optimal, :Suboptimal, :UserLimit]
    status_reroute = [:Infeasible]

    if m.status[:local_solve] in status_pass
        m.loglevel > 0 && println("local solver returns feasible point")
        bound_tightening(m, use_bound = true)    # performs bound-tightening with the local solve objective value
        (m.presolve_bt) && init_disc(m)          # Reinitialize discretization dictionary on tight bounds
        (m.disc_ratio_branch) && (m.disc_ratio = update_disc_ratio(m, true))
        add_partition(m, use_solution=m.best_sol)  # Setting up the initial discretization
        m.loglevel > 0 && println("presolve ended.")
    elseif m.status[:local_solve] in status_reroute
        (m.loglevel > 0) && println("first attempt at local solve failed, performing bound tightening without objective value...")
        bound_tightening(m, use_bound = false)                      # do bound tightening without objective value
        (m.loglevel > 0) && println("second attempt at local solve using tightened bounds...")
        local_solve(m, presolve = true) # local_solve(m) to generate a feasible solution which is a starting point for bounding_solve
        if m.status in status_pass  # successful second try
            (m.disc_ratio_branch) && (m.disc_ratio = update_disc_ratio(m, true))
            add_partition(m, use_solution=m.best_sol)
        else    # if this does not produce an feasible solution then solve atmc without discretization and use as a starting point
            (m.loglevel > 0) && println("reattempt at local solve failed, initialize discretization with lower bound solution... \n local solve remains infeasible...")
            create_bounding_mip(m)       # Build the bounding ATMC model
            bounding_solve(m)            # Solve bounding model
            (m.disc_ratio_branch) && (m.disc_ratio = update_disc_ratio(m))
            add_partition(m, use_solution=m.best_bound_sol)
        end
        m.loglevel > 0 && println("Presolve ended.")
    elseif m.status[:local_solve] == :Not_Enough_Degrees_Of_Freedom
        warn("Presolve ends with local solver yielding $(m.status[:local_solve]). \n Consider more replace equality constraints with >= and <= to resolve this.")
    else
        warn("Presolve ends with local solver yielding $(m.status[:local_solve]).")
    end

    cputime_presolve = time() - start_presolve
    m.logs[:presolve_time] += cputime_presolve
    m.logs[:total_time] = m.logs[:presolve_time]
    m.logs[:time_left] -= m.logs[:presolve_time]
    (m.loglevel > 0) && println("Presolve time = $(@compat round.(m.logs[:total_time],2))s")

    return
end

"""
    A wrapper function that collects some automated solver adjustments within the main while loop.
"""
function algorithm_automation(m::PODNonlinearModel)

    m.disc_var_pick == 3 && update_disc_cont_var(m)
    m.int_cumulative_disc && update_disc_int_var(m)

    if m.disc_ratio_branch
        m.disc_ratio = update_disc_ratio(m)    # Only perform for a maximum three times
    end

    return
end

"""
    Summarized function to determine whether to interrupt the main while loop.
"""
function check_exit(m::PODNonlinearModel)

    # Infeasibility check
    m.status[:bounding_solve] == :Infeasible && return true

    # Unbounded check
    m.status[:bounding_solve] == :Unbounded && return true

    # Optimality check
    m.best_rel_gap <= m.relgap && return true
    m.logs[:n_iter] >= m.maxiter && return true
    m.best_abs_gap <= m.absgap && return true

    # Userlimits check
    m.logs[:time_left] < m.tol && return true

    return false
end

"""

    local_solve(m::PODNonlinearModel, presolve::Bool=false)

Perform a local NLP or MINLP solve to obtain a feasible solution.
The `presolve` option is set to `true` when the function is invoked in [`presolve`](@ref).
Otherwise, the function is invoked from [`bounding_solve`](@ref).

"""

function local_solve(m::PODNonlinearModel; presolve = false)

    convertor = Dict(:Max=>:>, :Min=>:<)
    local_nlp_status = :Unknown
    do_heuristic = false

    var_type_screener = [i for i in m.var_type_orig if i in [:Bin, :Int]]

    if presolve
        if !isempty(var_type_screener) && m.minlp_solver != UnsetSolver()
            local_solve_model = interface_init_nonlinear_model(m.minlp_solver)
        elseif !isempty(var_type_screener)
            local_solve_model = interface_init_nonlinear_model(m.nlp_solver)
        else
            local_solve_model = interface_init_nonlinear_model(m.nlp_solver)
        end
    else
        local_solve_model = interface_init_nonlinear_model(m.nlp_solver)
    end

    if presolve == false
        l_var, u_var = fix_domains(m)
    else
        # l_var, u_var = m.l_var_orig, m.u_var_orig  # change to l_var_tight/u_var_tight ?
        l_var, u_var = m.l_var_tight[1:m.num_var_orig], m.u_var_tight[1:m.num_var_orig]
    end

    start_local_solve = time()
    interface_load_nonlinear_model(m, local_solve_model, l_var, u_var)
    (!m.d_orig.want_hess) && interface_init_nonlinear_data(m)

    if presolve == false
        interface_set_warmstart(local_solve_model, m.best_sol[1:m.num_var_orig])
    else
        initial_warmval = []
        for i in 1:m.num_var_orig
            isnan(m.d_orig.m.colVal[i]) ? push!(initial_warmval, 0.0) : push!(initial_warmval, m.d_orig.m.colVal[i])
        end
        interface_set_warmstart(local_solve_model, initial_warmval)
    end

    # The only case when MINLP solver is actually used
    if presolve && !isempty(var_type_screener)
        m.minlp_solver == UnsetSolver() || interface_set_vartype(local_solve_model, m.var_type_orig)
        interface_optimize(local_solve_model)
        if m.minlp_solver == UnsetSolver()
            do_heuristic = true
            local_nlp_status = :Heuristics
        else
            local_nlp_status = interface_get_status(local_solve_model)
            do_heuristic = false
        end
    else
        interface_optimize(local_solve_model)
        local_nlp_status = interface_get_status(local_solve_model)
    end

    cputime_local_solve = time() - start_local_solve
    m.logs[:total_time] += cputime_local_solve
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

    status_pass = [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
    status_heuristic = [:Heuristics]
    status_reroute = [:Infeasible]

    if local_nlp_status in status_pass
        candidate_obj = interface_get_objval(local_solve_model)
        candidate_sol = round.(interface_get_solution(local_solve_model), 5)
        update_incumb_objective(m, candidate_obj, candidate_sol)
        m.status[:local_solve] = local_nlp_status
        return
    elseif local_nlp_status in status_heuristic && do_heuristic
        m.status[:local_solve] = heu_basic_rounding(m, local_solve_model)
        return
    elseif local_nlp_status == :Infeasible
        heu_pool_multistart(m) == :LocalOptimal && return
        push!(m.logs[:obj], "INF")
        m.status[:local_solve] = :Infeasible
        return
    elseif local_nlp_status == :Unbounded
        push!(m.logs[:obj], "U")
        m.status[:local_solve] = :Unbounded
        presolve == true ? warn("[PRESOLVE] NLP local solve is unbounded.") : warn("[LOCAL SOLVE] NLP local solve is unbounded.")
        return
    else
        push!(m.logs[:obj], "E")
        m.status[:local_solve] = :Error
		presolve == true ? warn("[PRESOLVE] NLP solve failure $(local_nlp_status).") : warn("[LOCAL SOLVE] NLP local solve failure.")
        return
    end

    return
end


"""

    bounding_solve(m::PODNonlinearModel; kwargs...)

This is a solving process usually deal with a MIP or MIQCP problem for lower bounds of problems.
It solves the problem built upon a convexification base on a discretization Dictionary of some variables.
The convexification utilized is Tighten McCormick scheme.
See `create_bounding_mip` for more details of the problem solved here.

"""
function bounding_solve(m::PODNonlinearModel)

    convertor = Dict(:Max=>:<, :Min=>:>)
    boundlocator = Dict(:Max=>:+, :Min=>:-)
    boundlocator_rev = Dict(:Max=>:-, :Max=>:+)

    # Updates time metric and position solver
    update_mip_time_limit(m)
    update_boundstop_options(m)

    # ================= Solve Start ================ #
    start_bounding_solve = time()
    status = solve(m.model_mip, suppress_warnings=true)
    m.logs[:total_time] += time() - start_bounding_solve
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])
    # ================= Solve End ================ #

    status_solved = [:Optimal, :UserObjLimit, :UserLimit, :Suboptimal]
    status_maynosolution = [:UserObjLimit, :UserLimit]  # Watch out for these cases
    status_infeasible = [:Infeasible, :InfeasibleOrUnbounded]
    if status in status_solved
        (status == :Optimal) ? candidate_bound = m.model_mip.objVal : candidate_bound = m.model_mip.objBound
        candidate_bound_sol = [round.(getvalue(Variable(m.model_mip, i)), 6) for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
        if m.disc_consecutive_forbid > 0
            m.bound_sol_history[mod(m.logs[:n_iter]-1, m.disc_consecutive_forbid)+1] = copy(candidate_bound_sol) # Requires proper offseting
        end
        push!(m.logs[:bound], candidate_bound)
        if eval(convertor[m.sense_orig])(candidate_bound, m.best_bound)
            m.best_bound = candidate_bound
            m.best_bound_sol = copy(candidate_bound_sol)
            m.status[:bounding_solve] = status
            m.status[:bound] = :Detected
        end
        collect_lb_pool(m)    # Always collect details sub-optimal solution
    elseif status in status_infeasible
        push!(m.logs[:bound], "-")
        m.status[:bounding_solve] = :Infeasible
        PODDEBUG && print_iis_gurobi(m.model_mip) # Diagnostic code
        warn("[INFEASIBLE] Infeasibility detected via convex relaxation Infeasibility")
    elseif status == :Unbounded
        m.status[:bounding_solve] = :Unbounded
        warn("[UNBOUNDED] MIP solver return unbounded")
    else
        error("[EXCEPTION] MIP solver failure $(status)")
    end

    return
end

"""
    pick_disc_vars(m::PODNonlinearModel)

This function helps pick the variables for discretization. The method chosen depends on user-inputs.
In case when `indices::Int` is provided, the method is chosen as built-in method. Currently,
there exist two built-in method:

    * `max-cover(m.disc_var_pick=0, default)`: pick all variables involved in the non-linear term for discretization
    * `min-vertex-cover(m.disc_var_pick=1)`: pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified

For advance usage, `m.disc_var_pick` allows `::Function` inputs. User is required to perform flexible methods in choosing the non-linear variable.
For more information, read more details at [Hacking Solver](@ref).

"""
function pick_disc_vars(m::PODNonlinearModel)

    if isa(m.disc_var_pick, Function)
        eval(m.disc_var_pick)(m)
        length(m.disc_vars) == 0 && length(m.nonconvex_terms) > 0 && error("[USER FUNCTION] must select at least one variable to perform discretization for convexificiation purpose")
    elseif isa(m.disc_var_pick, Int)
        if m.disc_var_pick == 0
            ncvar_collect_nodes(m)
        elseif m.disc_var_pick == 1
            min_vertex_cover(m)
        elseif m.disc_var_pick == 2
            (length(m.candidate_disc_vars) > 15) ? min_vertex_cover(m) : ncvar_collect_nodes(m)
        elseif m.disc_var_pick == 3 # Initial
            (length(m.candidate_disc_vars) > 15) ? min_vertex_cover(m) : ncvar_collect_nodes(m)
        else
            error("Unsupported default indicator for picking variables for discretization")
        end
    else
        error("Input for parameter :disc_var_pick is illegal. Should be either a Int for default methods indexes or functional inputs.")
    end

    return
end
