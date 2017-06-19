"""
    Index function of the bound tightening method
"""
function presolve_bounds_tightening(m::PODNonlinearModel; kwargs...)

    if !m.presolve_do_bound_tightening
        return
    end

    if m.presolve_bound_tightening_method == 1
        minmax_bounds_tightening(m)         # Conduct basic min-max bound tightening scheme
    elseif m.presolve_bound_tightening_method == 2
        minmax_bounds_tightening(m, use_tmc=true)
    elseif isa(m.presolve_bound_tightening_method, Function)
        eval(m.presolve_bound_tightening_method)(m)
    else
        error("Unrecognized bound tightening scheme")
    end

    return
end

"""
    Sequential version of the bound tightening procedure.
    Algorithm layout goes here...
    Main implementation target...
"""
function minmax_bounds_tightening(m::PODNonlinearModel; kwargs...)

    # Some functinal constants
    both_senses = [:Min, :Max]             # Senses during bound tightening procedures
    tell_side = Dict(:Min=>1, :Max=>2)     # Positional information

    options = Dict(kwargs)

    # Regulating Speical Input Conditions: default use best feasible solution objective value
    haskey(options, :use_bound) ? bound = options[:use_bound] : bound = m.best_obj
    discretization = to_discretization(m, m.l_var_tight, m.u_var_tight)
    haskey(options, :use_tmc) ? discretization = add_discretization(m, use_solution=m.best_sol, use_discretization=discretization) : discretization = discretization
    discretization = resolve_lifted_var_bounds(m.nonlinear_info, discretization)
    # TODO: Potential risk above :: TMC with no feasible solution

    (m.log_level > 0) && println("start tightening bounds ...")
    (m.log_level > 99) && [println("[DEBUG] VAR $(var_idx) Original Bound [$(round(m.l_var_tight[var_idx],4)) < - > $(round(m.u_var_tight[var_idx],4))]") for var_idx in m.all_nonlinear_vars]

    # ======= Algorithm Starts ======= #
    keeptighening = true
    while keeptighening && m.logs[:time_left] > m.tolerance && m.logs[:bt_iter] <= m.presolve_maxiter # Stopping criteria

        keeptighening = false
        m.logs[:bt_iter] += 1
        (m.log_level > 99) && println("[DEBUG] Iteration - $(m.logs[:bt_iter])")
        temp_bounds = Dict()

        # Perform Bound Contraction
        for var_idx in m.all_nonlinear_vars
            temp_bounds[var_idx] = [discretization[var_idx][1], discretization[var_idx][end]]
            create_bound_tightening_model(m, discretization, bound)
            for sense in both_senses
                @objective(m.model_mip, sense, Variable(m.model_mip, var_idx))
                solve_bound_tightening_model(m)
                temp_bounds[var_idx][tell_side[sense]] = getobjectivevalue(m.model_mip)
                m.log_level > 99 && println("[DEBUG] contracting VAR $(var_idx) with $(sense) problem, results in $(getobjectivevalue(m.model_mip)) from $(temp_bounds[var_idx])")
            end
        end

        # Updates the discretizatio structure
        for var_idx in m.all_nonlinear_vars
            if abs((temp_bounds[var_idx][1] - discretization[var_idx][1])/discretization[var_idx][1]) > m.presolve_tolerance
                keeptighening = true # Continue to perform the next iteration
                discretization[var_idx][1] = temp_bounds[var_idx][1]
            end
            if abs((discretization[var_idx][end]-temp_bounds[var_idx][end])/discretization[var_idx][end]) > m.presolve_tolerance
                (m.log_level > 0) && print("+")
                keeptighening = true
                discretization[var_idx][end] = temp_bounds[var_idx][2]
            end
        end

        discretization = resolve_lifted_var_bounds(m.nonlinear_info, discretization)
        haskey(options, :use_tmc) ? discretization = add_discretization(m, use_solution=m.best_sol, use_discretization=flatten_discretization(discretization)) : discretization = discretization
    end

    (m.log_level > 0) && println("\nfinished bound tightening in $(m.logs[:bt_iter])iterations, applying tighten bounds")
    # ======= Algorithm Ends ======== #

    # Updae the bounds with the tighten ones
    m.l_var_tight, m.u_var_tight = update_var_bounds(discretization)

    (m.log_level > 99)  && [println("[DEBUG] VAR $(i) BOUND contracted |$(round(m.l_var_orig[i],4)) --> | $(round(m.l_var_tight[i],4)) - * - $(round(m.u_var_tight[i],4)) | <-- $(round(m.u_var_orig[i],4)) |") for i in m.all_nonlinear_vars]
    (m.log_level > 0) && print("\n")
    return
end

"""
    This function takes in the initial discretization information and builds a bound tighting model that is connected to .model_mip
    It is an algorithm specific function
"""
function create_bound_tightening_model(m::PODNonlinearModel, discretization, bound; kwargs...)

    options = Dict(kwargs)

    start_build = time()

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP model
    post_amp_vars(m, use_discretization=discretization)
    post_amp_lifted_constraints(m)
    post_amp_mccormick(m, use_discretization=discretization)
    post_obj_bounds(m, bound)

    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build * m.presolve_track_time
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time] * m.presolve_track_time)

    return
end

function solve_bound_tightening_model(m::PODNonlinearModel; kwargs...)

    # ========= MIP Solve ========= #
    if m.presolve_mip_timelimit < Inf
        update_mip_time_limit(m, timelimit = max(0.0, min(m.presolve_mip_timelimit, m.timeout - m.logs[:total_time])))
    else
        update_mip_time_limit(m, timelimit = max(0.0, m.timeout - m.logs[:total_time]))
    end

    start_solve = time()
    status = solve(m.model_mip, suppress_warnings=true, relaxation=m.presolve_mip_relaxation)
    cputime_solve = time() - start_solve
    m.logs[:total_time] += cputime_solve * m.presolve_track_time
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time] * m.presolve_track_time)
    # ========= MIP Solve ========= #

    return
end


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

function resolve_closed_var_bounds(m::PODNonlinearModel; kwargs...)

    for var in m.all_nonlinear_vars
        if abs(m.l_var_tight[var] - m.u_var_tight[var]) < m.tolerance   # Closed Bound Criteria
            # Clean nonlinear_info by deleting the info

        end
    end

    return
end
