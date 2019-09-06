"""
    High-level Function
"""
function MathProgBase.optimize!(m::AlpineNonlinearModel)
    if m.presolve_infeasible
        summary_status(m)
        return
    end
    presolve(m)
    global_solve(m)
    m.loglevel > 0 && logging_row_entry(m, finsih_entry=true)
    println("====================================================================================================")
    summary_status(m)  
    return
end

"""
    global_solve(m::AlpineNonlinearModel)

Perform global optimization algorithm that is based on the adaptive piecewise convexification.
This iterative algorithm loops over [`bounding_solve`](@ref) and [`local_solve`](@ref) until the optimality gap between the lower bound (relaxed problem with min. objective) and the upper bound (feasible problem) is within the user prescribed limits.
Each [`bounding_solve`](@ref) provides a lower bound that serves as the partioning point for the next iteration (this feature can be modified given a different `add_adaptive_partition`).
Each [`local_solve`](@ref) provides an incumbent feasible solution. The algorithm terminates when atleast one of these conditions are satisfied: time limit, optimality condition, or iteration limit.

"""
function global_solve(m::AlpineNonlinearModel)

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

function run_bounding_iteration(m::AlpineNonlinearModel)
    m.loglevel > 0 && logging_head(m)
    m.presolve_track_time || reset_timer(m)

    m.logs[:n_iter] += 1
    create_bounding_mip(m)                  # Build the relaxation model
    bounding_solve(m)                       # Solve the relaxation model
    update_opt_gap(m)                       # Update optimality gap
    check_exit(m) && return                 # Feasibility check
    m.loglevel > 0 && logging_row_entry(m)  # Logging
    local_solve(m)                          # Solve local model for feasible solution
    update_opt_gap(m)                       # Update optimality gap
    check_exit(m) && return                 # Detect optimality termination
    algorithm_automation(m)                 # Automated adjustments
    add_partition(m)                        # Add extra discretizations


    return
end 

"""
    presolve(m::AlpineNonlinearModel)
"""
function presolve(m::AlpineNonlinearModel)

    start_presolve = time()
    m.loglevel > 0 && printstyled("PRESOLVE \n", color=:cyan)
    m.loglevel > 0 && println("  Doing local search")
    local_solve(m, presolve = true)

    # Solver status - returns error when see different
    status_pass = [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
    status_reroute = [:Infeasible, :Infeasibles]

    if m.status[:local_solve] in status_pass
        m.loglevel > 0 && println("  Local solver returns a feasible point")
        bound_tightening(m, use_bound = true)    # performs bound-tightening with the local solve objective value
        m.presolve_bt && init_disc(m)            # Re-initialize discretization dictionary on tight bounds
        m.disc_ratio_branch && (m.disc_ratio = update_disc_ratio(m, true))
        add_partition(m, use_solution=m.best_sol)  # Setting up the initial discretization
        # m.loglevel > 0 && println("Ending the presolve")
    elseif m.status[:local_solve] in status_reroute
        (m.loglevel > 0) && println("  Bound tightening without objective bounds (OBBT)")
        bound_tightening(m, use_bound = false)                      # do bound tightening without objective value
        (m.disc_ratio_branch) && (m.disc_ratio = update_disc_ratio(m))
        m.presolve_bt && init_disc(m)
        # m.loglevel > 0 && println("Ending the presolve")
    elseif m.status[:local_solve] == :Not_Enough_Degrees_Of_Freedom
        @warn " Warning: Presolve ends with local solver yielding $(m.status[:local_solve]). \n Consider more replace equality constraints with >= and <= to resolve this."
    else
        @warn " Warning: Presolve ends with local solver yielding $(m.status[:local_solve])."
    end

    cputime_presolve = time() - start_presolve
    m.logs[:presolve_time] += cputime_presolve
    m.logs[:total_time] = m.logs[:presolve_time]
    m.logs[:time_left] -= m.logs[:presolve_time]
    # (m.loglevel > 0) && println("Presolve time = $(round.(m.logs[:total_time]; digits=2))s")
    (m.loglevel > 0) && println("  Completed presolve in $(round.(m.logs[:total_time]; digits=2))s ($(m.logs[:bt_iter]) iterations).")
    return
end

"""
    A wrapper function that collects some automated solver adjustments within the main while loop.
"""
function algorithm_automation(m::AlpineNonlinearModel)

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
function check_exit(m::AlpineNonlinearModel)

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
    local_solve(m::AlpineNonlinearModel, presolve::Bool=false)

Perform a local NLP or MINLP solve to obtain a feasible solution.
The `presolve` option is set to `true` when the function is invoked in [`presolve`](@ref).
Otherwise, the function is invoked from [`bounding_solve`](@ref).
"""
function local_solve(m::AlpineNonlinearModel; presolve = false)

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
        l_var, u_var = m.l_var_tight[1:m.num_var_orig], m.u_var_tight[1:m.num_var_orig]
    end

    start_local_solve = time()
    interface_load_nonlinear_model(m, local_solve_model, l_var, u_var)
    (!m.d_orig.want_hess) && interface_init_nonlinear_data(m.d_orig)

    if presolve == false
        interface_set_warmstart(local_solve_model, m.best_sol[1:m.num_var_orig])
    else
        initial_warmval = Float64[]
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
    status_reroute = [:Infeasible, :Infeasibles]
    
    if local_nlp_status in status_pass
        candidate_obj = interface_get_objval(local_solve_model)
        candidate_sol = round.(interface_get_solution(local_solve_model); digits=5)
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
        if presolve == true 
            @warn "  Warning: NLP local solve is unbounded." 
        else 
            @warn "  Warning: NLP local solve is unbounded." 
        end
        return
    else
        push!(m.logs[:obj], "E")
        m.status[:local_solve] = :Error
        if presolve == true 
            @warn " Warning: NLP solve failure $(local_nlp_status)." 
        else 
            @warn " Warning: NLP local solve failure."
        end
        return
    end

    return
end


"""

    bounding_solve(m::AlpineNonlinearModel; kwargs...)

This is a solving process usually deal with a MIP or MIQCP problem for lower bounds of problems.
It solves the problem built upon a convexification base on a discretization Dictionary of some variables.
The convexification utilized is Tighten McCormick scheme.
See `create_bounding_mip` for more details of the problem solved here.

"""
function bounding_solve(m::AlpineNonlinearModel)

    convertor = Dict(:Max=>:<, :Min=>:>)
    boundlocator = Dict(:Max=>:+, :Min=>:-)
    boundlocator_rev = Dict(:Max=>:-, :Max=>:+)

    # Updates time metric and the termination bounds
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
        candidate_bound_sol = [round.(getvalue(Variable(m.model_mip, i)); digits=6) for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
        # Experimental code
        measure_relaxed_deviation(m, sol=candidate_bound_sol)
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
        ALPINE_DEBUG && print_iis_gurobi(m.model_mip) # Diagnostic code
        @warn "  Warning: Infeasibility detected via convex relaxation Infeasibility"
    elseif status == :Unbounded
        m.status[:bounding_solve] = :Unbounded
        @warn "  Warning: MIP solver return unbounded"
    else
        error("  Warning: MIP solver failure $(status)")
    end

    return
end

"""
    pick_disc_vars(m::AlpineNonlinearModel)

This function helps pick the variables for discretization. The method chosen depends on user-inputs.
In case when `indices::Int` is provided, the method is chosen as built-in method. Currently,
there are two built-in options for users as follows:

    * `max-cover (m.disc_var_pick=0, default)`: pick all variables involved in the non-linear term for discretization
    * `min-vertex-cover (m.disc_var_pick=1)`: pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified

For advanced usage, `m.disc_var_pick` allows `::Function` inputs. User can provide his/her own function to choose the variables for discretization.

"""
function pick_disc_vars(m::AlpineNonlinearModel)

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
