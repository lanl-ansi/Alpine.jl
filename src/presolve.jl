"""
    Index function of the bound tightening method
"""
function presolve_bounds_tightening(m::PODNonlinearModel; kwargs...)

    if !m.do_bound_tightening
        return
    end

    if m.bound_tightening_method == 1
        minmax_bounds_tightening(m)         # Conduct basic min-max bound tightening scheme
    elseif m.bound_tightening_method == 2
        minmax_bounds_tightening(m, use_tmc=true)
    elseif isa(m.bound_tightening_method, Function)
        eval(m.bound_tightening_method)(m)
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
    bt_discretization = copy_discretization(m, m.l_var_tight, m.u_var_tight)
    haskey(options, :use_tmc) ? bt_discretization = add_discretization(m, use_solution=m.best_sol, use_discretization=bt_discretization) : bt_discretization = bt_discretization
    # Potential risk above :: TMC with no feasible solution

    @show bt_discretization
    @show m.bound_tightening_method

    m.log_level > 0 &&
        println("start tightening bounds ...")

    # ======= Algorithm Starts ======= #
    keeptighening = true
    while keeptighening && m.logs[:time_left] > m.tolerance # Stopping criteria
        keeptighening = false
        (m.log_level > 99) &&
            println("Iteration -$(m.logs[:bt_iter])")

        for var_idx in m.all_nonlinear_vars
            (m.log_level > 99 &&
                m.logs[:bt_iter]==0) && println("[VAR $(var_idx)] Original Bound [$(round(m.l_var_tight[var_idx],4)) < - > $(round(m.u_var_tight[var_idx],4))]")

            temp_bounds = [bt_discretization[var_idx][1], bt_discretization[var_idx][end]]
            create_bound_tightening_model(m, bt_discretization, bound)
            for sense in both_senses
                @objective(m.model_mip, sense, Variable(m.model_mip, var_idx))
                status = solve(m.model_mip, suppress_warnings=true)
                temp_bounds[tell_side[sense]] = getobjectivevalue(m.model_mip)
            end

            print(m.model_mip)
            error("Stop")
            m.log_level > 99 &&
                println("[VAR $(var_idx)] Bound [$(round(temp_bounds[1],4)) < - > $(round(temp_bounds[end],4))]")

            if abs((temp_bounds[1] - bt_discretization[var_idx][1])/bt_discretization[var_idx][1]) > 1e-3 || abs((bt_discretization[var_idx][end]-temp_bounds[end])/bt_discretization[var_idx][end] > 1e-3)
                (m.log_level > 0) && print("+")
                keeptighening = true # Continue to perform the next iteration
                bt_discretization[var_idx][1] = max(temp_bounds[1], bt_discretization[var_idx][1] + m.tolerance)
                bt_discretization[var_idx][end] = min(temp_bounds[2], bt_discretization[var_idx][end] - m.tolerance)
                bt_discretization = resolve_lifted_var_bounds(m.nonlinear_info, bt_discretization)
            end
        end
        (m.log_level > 99) && println(" -> ")
        m.logs[:bt_iter] += 1
    end

    (m.log_level > 0) && println("applying tighten bounds")
    m.l_var_tight, m.u_var_tight = update_var_bounds(bt_discretization)
    # ======= Algorithm Ends ======== #

    (m.log_level > 0) && print("\n")
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

"""
    This function takes in the initial discretization information and builds a bound tighting model that is connected to .model_mip
    It is an algorithm specific function
"""
function create_bound_tightening_model(m::PODNonlinearModel, discretization, bound; kwargs...)

    options = Dict(kwargs)

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP model
    post_amp_vars(m, use_discretization=discretization)
    post_amp_lifted_constraints(m)
    post_amp_mccormick(m, use_discretization=discretization)
    post_obj_bounds(m, bound)

    return
end

"""
    Take the objective expression and post an limitation bound on it.
"""
function post_obj_bounds(m::PODNonlinearModel, bound::Float64; kwargs...)
    if m.sense_orig == :Max
        @constraint(m.model_mip,
            sum(m.lifted_obj_aff_mip[:coefs][j]*Variable(m.model_mip, m.lifted_obj_aff_mip[:vars][j].args[2]) for j in 1:m.lifted_obj_aff_mip[:cnt]) >= bound)
    elseif m.sense_orig == :Min
        @constraint(m.model_mip,
            sum(m.lifted_obj_aff_mip[:coefs][j]*Variable(m.model_mip, m.lifted_obj_aff_mip[:vars][j].args[2]) for j in 1:m.lifted_obj_aff_mip[:cnt]) <= bound)
    end

    return
end
