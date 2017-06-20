"""

    presolve_bounds_tightening(m::PODNonlinearModel)

High-level presolve caller that activate the presolve process. The built-in presolve
procedure is bound tightening method in constraint programming. The goal of these built-in
methods is to contract variable feasible domains using constraint programming scheme.

Currently, two bounding tightening method is implemented [`minmax_bounds_tightening`](@ref).

    * Min-Max McCormick Bounds Tightening
    * Min-Max Tighten McCormick Bounds Tightening

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

    minmax_bounds_tightening(m:PODNonlinearModel; use_bound::Float64, use_tmc::Bool)

This function implements the algorithm used to tight variable feasible domains using constraint programming shceme.
It utilize McCormick or Tighten McCormick Realxation to creat MILP and solve by optimizing on variable bounds.

The algorithm iteratively walk through all variables involved in non-linear terms (discretizating variables) and solve
min-max problems for optimal lower-upper bounds. During the procedure, any tighten bounds will be utilized for further
bound tightening operations. For better bounds contraction quality, `use_tmc` can be assigned `true` to use a slightly
tighter formulation while in sacrificing the solution time. The bound tightening algirthm terminates base on improvements,
time limits, and iteration limits.

It is highly recommended that `use_bound` is available for best quality ensurance. Several other options is available
for the algorithm tuning. For more details, see [Parameters](@ref).

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
    create_bound_tightening_model(m::PODNonlinearModel, discretization::Dict, bound::Float64)

This function takes in the initial discretization information and builds a bound tighting model that is connected to .model_mip
It is an algorithm specific function in [`minmax_bounds_tightening`](@ref) for best felxibility in tuning and hacking.

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

"""

    solve_bound_tightening_model(m::PODNonlinearModels)

A function that solves the min-max model used in built-in presolve algorithm.

"""
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

"""
    resolve_lifted_var_bounds(nonlinear_info::Dict, discretization::Dict)

For discretization to be performed, it is not allowed for a discretizing variable to have inifinate bounds.
This violation can be seen in new lifted proposed for multi-variant non-linear terms. This functions resolve
these issues by use the problem information to reason some lifted variables bounds. Such bounds reasoning
can also speed up the bounding problem during the main algorithm by providing trivial information to the solver.
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

function resolve_closed_var_bounds(m::PODNonlinearModel; kwargs...)

    for var in m.all_nonlinear_vars
        if abs(m.l_var_tight[var] - m.u_var_tight[var]) < m.tolerance   # Closed Bound Criteria
            # Clean nonlinear_info by deleting the info

        end
    end

    return
end
