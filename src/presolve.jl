"""
    bound_tightening(m::Optimizer)

Entry point to the optimization-based bound-tightening (OBBT) algorithm. The aim of the OBBT algorithm
is to sequentially tighten the variable bounds until a fixed point is reached.

Currently, two OBBT methods are implemented in [`minmax_bound_tightening`](@ref).

    * Bound-tightening with polyhedral relaxations (McCormick, Lambda for convex-hull)
    * Bound-tightening with piecewise polyhedral relaxations: (with three partitions around the local feasible solution)
If no local feasible solution is obtained, the algorithm defaults to OBBT without partitions
"""
function bound_tightening(m::Optimizer; use_bound = true, kwargs...)

    Alp.get_option(m, :presolve_bt) || return

    if Alp.get_option(m, :presolve_bt_algo) == 1
        Alp.minmax_bound_tightening(m, use_bound=use_bound)
    elseif Alp.get_option(m, :presolve_bt_algo) == 2
        Alp.minmax_bound_tightening(m, use_bound=use_bound, use_tmc=true)
    elseif isa(Alp.get_option(m, :presolve_bt_algo), Function)
        # eval(Alp.get_option(m, :presolve_bt_algo))(m)
        Alp.get_option(m, :presolve_bt_algo)(m)
    else
        error("Unrecognized bound tightening algorithm")
    end

    return
end

"""
    minmax_bound_tightening(m:Optimizer; use_bound::Bool=true, use_tmc::Bool)

This function implements the OBBT algorithm to tighten the variable bounds.
It utilizes either the basic polyhedral relaxations or the piecewise polyhedral relaxations (TMC)
to tighten the bounds. The TMC has additional binary variables while performing OBBT.

The algorithm as two main parameters. The first is the `use_tmc`, which when set to `true`
invokes the algorithm on the TMC relaxation. The second parameter `use_bound` takes in the
objective value of the local solve solution stored in `best_sol` for performing OBBT. The `use_bound` option is set
to `true` when the local solve is successful in obtaining a feasible solution, else this parameter
is set to `false`.

For details, refer to section 3.1.1 of
Nagarajan, Lu, Wang, Bent, Sundar, "An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs"
[link](https://doi.org/10.1007/s10898-018-00734-1).

Several other user-input options can be used to tune the OBBT algorithm.
For more details, see [Presolve Options](https://lanl-ansi.github.io/Alpine.jl/latest/parameters/#Presolve-Options).

"""
function minmax_bound_tightening(m::Optimizer; use_bound = true, time_limit = Inf, kwargs...)

    # Some functinal constants
    both_senses = [MOI.MIN_SENSE, MOI.MAX_SENSE]             # Senses during bound tightening procedures
    tell_side = Dict(MOI.MIN_SENSE => 1, MOI.MAX_SENSE=>2)     # Positional information
    tell_round = Dict(MOI.MIN_SENSE => floor, MOI.MAX_SENSE => ceil)

    options = Dict(kwargs)

    st = time() # Track start time
    if time_limit == Inf
        time_limit = Alp.get_option(m, :presolve_bt_time_limit)
    end

    # Regulating special input conditions: default use best feasible solution objective value
    (use_bound == true) ? bound = m.best_obj : bound = Inf
    l_var_orig = copy(m.l_var_tight)
    u_var_orig = copy(m.u_var_tight)

    discretization = Alp.to_discretization(m, m.l_var_tight, m.u_var_tight)
    if use_bound == false && haskey(options, :use_tmc)
        (Alp.get_option(m, :log_level) > 0) && @warn " Local solve infeasible; defaulting to doing bound-tightening without partitions."
    end
    if use_bound == true && haskey(options, :use_tmc)
        discretization = Alp.add_adaptive_partition(m, use_solution=m.best_sol, use_disc=discretization)
    end
    discretization = Alp.resolve_var_bounds(m, discretization) # recomputation of bounds for lifted_variables

    (Alp.get_option(m, :log_level) > 0) && println("  Starting bound-tightening")

    width_tol  = Alp.get_option(m, :presolve_bt_width_tol)
    bound_tol  = Alp.get_option(m, :presolve_bt_bound_tol)
    improv_tol = Alp.get_option(m, :presolve_bt_improv_tol)

    keep_tightening = true

    # start of the solve
    while keep_tightening && (m.logs[:time_left] > Alp.get_option(m, :tol)) && (m.logs[:bt_iter] < Alp.get_option(m, :presolve_bt_max_iter)) # Stopping criteria

        keep_tightening = false
        m.logs[:bt_iter] += 1
        Alp.get_option(m, :log_level) > 199 && println("  Iteration - $(m.logs[:bt_iter])")
        temp_bounds = Dict()

        avg_reduction   = 0.0
        max_reduction   = 0.0
        total_reduction = 0.0
        max_width       = 0.0
        current_rel_gap = 0.0

        # Perform sequential bound tightening
        for var_idx in m.candidate_disc_vars

            temp_bounds[var_idx] = [discretization[var_idx][1], discretization[var_idx][end]]
            if (discretization[var_idx][end] - discretization[var_idx][1]) > width_tol
                Alp.create_bound_tightening_model(m, discretization, bound)
                for sense in both_senses
                    JuMP.@objective(m.model_mip, sense, _index_to_variable_ref(m.model_mip, var_idx))
                    stats = Alp.solve_bound_tightening_model(m)
                    if stats["status"] in STATUS_OPT
                        temp_bounds[var_idx][tell_side[sense]] = tell_round[sense](JuMP.objective_value(m.model_mip)/bound_tol)*bound_tol  # Objective truncation for numerical issues
                    elseif stats["status"] in STATUS_LIMIT
                        temp_bounds[var_idx][tell_side[sense]] = tell_round[sense](JuMP.objective_bound(m.model_mip)/bound_tol)*bound_tol
                    elseif stats["status"] in STATUS_INF
                        @warn("Infeasible model detected within bound tightening - bounds not updated")
                    else
                        @warn("Unknown status within  bound tightening models")
                    end
                end
            end

            if (temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]] > temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]])
                temp_bounds[var_idx] = [discretization[var_idx][1], discretization[var_idx][end]]
            end
            if (temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]] > discretization[var_idx][end])
                temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]] = discretization[var_idx][1]
            end
            if (temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] < discretization[var_idx][1])
                temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] = discretization[var_idx][end]
            end

            bound_reduction = 0.0
            if (temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] - temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]]) <= width_tol
                midpoint = (temp_bounds[var_idx][1] + temp_bounds[var_idx][end]) / 2
                if (midpoint - discretization[var_idx][1] < width_tol/2)
                    temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]] = discretization[var_idx][1]
                    temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] = discretization[var_idx][1] + (width_tol)
                elseif (discretization[var_idx][end] - midpoint < width_tol/2)
                    temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]] = discretization[var_idx][end] - (width_tol)
                    temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] = discretization[var_idx][end]
                else
                    temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]] = midpoint - (width_tol/2)
                    temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] = midpoint + (width_tol/2)
                end
            end

            new_range = round(temp_bounds[var_idx][tell_side[MOI.MAX_SENSE]] - temp_bounds[var_idx][tell_side[MOI.MIN_SENSE]], digits = 4)
            old_range = discretization[var_idx][end] - discretization[var_idx][1]
            bound_reduction = old_range - new_range
            total_reduction += bound_reduction
            max_reduction = max(bound_reduction, max_reduction)
            max_width = max(new_range, max_width)

            (Alp.get_option(m, :log_level) > 99) && print("+")
            (Alp.get_option(m, :log_level) > 99) && println("  VAR $(var_idx) LB contracted $(discretization[var_idx][1])=>$(temp_bounds[var_idx][1])")
            (Alp.get_option(m, :log_level) > 99) && print("+")
            (Alp.get_option(m, :log_level) > 99) && println("  VAR $(var_idx) UB contracted $(discretization[var_idx][end])=>$(temp_bounds[var_idx][end])")
            
            discretization[var_idx][1]   = temp_bounds[var_idx][1]
            discretization[var_idx][end] = temp_bounds[var_idx][end]

        end

        avg_reduction = total_reduction/length(keys(temp_bounds))

        bound_avg_reduction = (avg_reduction > improv_tol)
        bound_max_reduction = (max_reduction > improv_tol)
        bound_max_width     = (max_width > width_tol)

        # Deactivate this termination criterion if it slows down the OBBT convergence 
        stats = Alp.relaxation_model_obbt(m, discretization, bound)
        if Alp.is_min_sense(m)
            current_rel_gap = Alp.eval_opt_gap(m, stats["relaxed_obj"], bound)
        elseif Alp.is_max_sense(m)
            current_rel_gap = Alp.eval_opt_gap(m, bound, stats["relaxed_obj"])
        end

        keep_tightening = (bound_avg_reduction) && (bound_max_reduction) && (bound_max_width) && (current_rel_gap > Alp.get_option(m, :rel_gap))

        if !keep_tightening && !isinf(current_rel_gap)
            m.presolve_best_rel_gap = current_rel_gap * 100
        end 

        discretization = Alp.resolve_var_bounds(m, discretization)
        
        if haskey(options, :use_tmc)
            discretization = Alp.add_adaptive_partition(m, use_solution=m.best_sol, use_disc = Alp.flatten_discretization(discretization))
        else
            discretization = discretization
        end

        time() - st > time_limit && break

    end

    (Alp.get_option(m, :log_level) > 1) && println("  Variables whose bounds were tightened:")
    (Alp.get_option(m, :log_level) > 0) && println("  Actual iterations (OBBT): ", (m.logs[:bt_iter]))
    
    m.l_var_tight, m.u_var_tight = Alp.update_var_bounds(discretization)
    m.discretization = Alp.add_adaptive_partition(m, use_solution = m.best_sol)

    for i in m.disc_vars
        contract_ratio = round(1-abs(m.l_var_tight[i] - m.u_var_tight[i])/abs(l_var_orig[i] - u_var_orig[i]); digits=2)*100
        if Alp.get_option(m, :log_level) > 0 && contract_ratio > 0.0001
            (Alp.get_option(m, :log_level) > 1) && (println("    VAR $(i): $(contract_ratio)% contraction |$(round(l_var_orig[i]; digits=4)) --> | $(round(m.l_var_tight[i]; digits=4)) - $(round(m.u_var_tight[i]; digits=4)) | <-- $(round(u_var_orig[i]; digits=4)) |"))
        end
    end
    # (Alp.get_option(m, :log_level) > 0) && print("\n")
    return
end

"""
    create_bound_tightening_model(m::Optimizer, discretization::Dict, bound::Float64)

This function takes in the initial discretization information and builds the OBBT model.
It is an algorithm specific function called by [`minmax_bound_tightening`](@ref).
"""
function create_bound_tightening_model(m::Optimizer, discretization, bound::Number; kwargs...)

    # options = Dict(kwargs)
    start_build = time()
    m.model_mip = Model(Alp.get_option(m, :mip_solver)) # Construct JuMP model
    Alp.amp_post_vars(m, use_disc = discretization)
    Alp.amp_post_lifted_constraints(m)
    Alp.amp_post_convexification(m, use_disc = discretization)  # Convexify problem

    if !(isinf(bound))
        Alp.post_objective_bound(m, bound)
    else 
        @warn "Dropping the objective bound constraint in presolve bound tightening"
    end

    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, Alp.get_option(m, :time_limit) - m.logs[:total_time])

    return
end

function relaxation_model_obbt(m::Optimizer, discretization, bound::Number)
    Alp.create_bound_tightening_model(m, discretization, bound)

    obj_expr = sum(m.bounding_obj_mip[:coefs][j]*_index_to_variable_ref(m.model_mip, m.bounding_obj_mip[:vars][j].args[2]) 
                   for j in 1:m.bounding_obj_mip[:cnt])

    if Alp.is_min_sense(m)
        JuMP.@objective(m.model_mip, Min, obj_expr)
    elseif Alp.is_max_sense(m)
        JuMP.@objective(m.model_mip, Max, obj_expr)
    end

    stats = Alp.solve_bound_tightening_model(m)
    
    return stats
end

"""

    solve_bound_tightening_model(m::Optimizer)

A function that solves the min and max OBBT model.

"""
function solve_bound_tightening_model(m::Optimizer; kwargs...)
    stats = Dict()

    # ========= MILP Solve ========= #
    time_limit = max(0.0, if Alp.get_option(m, :presolve_bt_mip_time_limit) < Inf
        min(Alp.get_option(m, :presolve_bt_mip_time_limit), Alp.get_option(m, :time_limit) - m.logs[:total_time])
    else
        Alp.get_option(m, :time_limit) - m.logs[:total_time]
    end)
    MOI.set(m.model_mip, MOI.TimeLimitSec(), time_limit)

    start_solve = time()
    if Alp.get_option(m, :presolve_bt_relax_integrality)
        JuMP.relax_integrality(m.model_mip)
    end

    JuMP.optimize!(m.model_mip)
    status = MOI.get(m.model_mip, MOI.TerminationStatus())

    stats["status"] = status
    stats["relaxed_obj"] = JuMP.objective_value(m.model_mip)
    
    cputime_solve = time() - start_solve
    m.logs[:total_time] += cputime_solve
    m.logs[:time_left] = max(0.0, Alp.get_option(m, :time_limit) - m.logs[:total_time])
    # ========= MILP Solve ========= #

    return stats
end

"""
    post_objective_bound(m::Optimizer, bound::Float64; kwargs...)

This function adds the upper/lower bounding constraint on the objective function 
for the optimization models solved within the OBBT algorithm. 
"""
function post_objective_bound(m::Optimizer, bound::Number; kwargs...)

    obj_expr = sum(m.bounding_obj_mip[:coefs][j]*_index_to_variable_ref(m.model_mip, m.bounding_obj_mip[:vars][j].args[2]) 
                   for j in 1:m.bounding_obj_mip[:cnt])
    
    obj_bound_tol = Alp.get_option(m, :presolve_bt_obj_bound_tol)

    if Alp.is_max_sense(m)
        JuMP.@constraint(m.model_mip, obj_expr >= floor(bound/obj_bound_tol) * obj_bound_tol)
    elseif Alp.is_min_sense(m)
        JuMP.@constraint(m.model_mip, obj_expr <=  ceil(bound/obj_bound_tol) * obj_bound_tol)
    end

    return
end
