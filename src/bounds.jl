@doc """
    Sequential version of the bound tightening procedure.
    Algorithm layout goes here...
    Main implementation target...
""" ->
function presolve_bounds_tightening(m::PODNonlinearModel; kwargs...)

    both_senses = [:Min, :Max]
    tell_bounds = Dict(:Min=>1, :Max=>2)     # Positional information

    # ======= Algorithm Starts ======= #
    keeptighening = true
    while keeptighening  # Stopping criteria
        new_bounds = Dict()
        # Build BT Model
        m.model_mip = Model(solver=m.mip_solver) # Construct jump model
        post_amp_vars(m)
        post_amp_lifted_constraints(m)  mn
        post_amp_mccormick(m, use_discretization=bt_discretization)

        for i = 1:length(m.all_nonlinear_vars)
            var_idx = m.all_nonlinear_vars[i]
            bt_discretization = m.discretization     # TODO
            new_bounds[var_idx] = [m.l_var_orig[var_idx], m.u_var_orig[var_idx]]
            for sense in both_senses
                # TODO : why do we need to do this multiple times?
                @objective(m.model_mip, sense, Variable(m.model_mip, var_idx))
                status = solve(m.model_mip, suppress_warnings=true)
                new_bounds[tell_bounds[sense]] = getobjectivevalue(m.model_mip)
            end
        end

        # Update ounds for lifted variables
        # for key in sort(collect(keys(NLmap)))  # update bounds for relaxed vars
        #    lb, ub = updatebounds(⋅)
        # end

        m.logs[:bt_iter] += 1
    end

    return
end


@doc """
    Initialize the discretization of all variables by setting up the bounds.
""" ->
function initialize_discretization(m::PODNonlinearModel;kwargs...)

    for var in 1:m.num_var_orig
        # TODO :: add a scheme for equal bounds
        lb = m.l_var_orig[var]
        ub = m.u_var_orig[var]
        if (lb == Inf) || (ub == -Inf)
            m.log_level > 99 && println("[VAR$(var)] Infinite bounds.")
        end
        m.discretization[var] = [lb, ub]
    end

    return
end

@doc """
    Basic method used to add a new partition.
    This method make modiciation in .discretization
    Keyward arguments:
        use_solution :: Vector{Float64}  -> Solution vector
""" ->
function add_discretization(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_solution) ? point_vec = options[:use_solution] : point_vec = m.best_bound_sol

    #===============================================================
    Consider original partition [0, 3, 7, 9], where LB/any solution is 4.
    Use ^ as the new partition, | as the original partition

    A case when discretize ratio = 4
    | -------- | - ^ -- * -- ^ ---- | -------- |
    0          3  3.5   4   4.5     7          9

    A special case when discretize ratio = 2
    | -------- | --- * --- ^ ---- | -------- |
    0          3     4     5      7          9
    ===============================================================#

    # ? Perform discretization base on type of nonlinear terms

    for i in 1:m.num_var_orig
        point = point_vec[i]
        if i in m.var_discretization_mip  # Only construct when discretized
            for j in 1:length(m.discretization[i])
                if point >= m.discretization[i][j] && point <= m.discretization[i][j+1]  # Locating the right location
                    @assert j < length(m.discretization[i])
                    lb_local = m.discretization[i][j]
                    ub_local = m.discretization[i][j+1]
                    distance = ub_local - lb_local
                    radius = distance / m.discretization_ratio
                    lb_new = max(point - radius, lb_local)
                    ub_new = min(point + radius, ub_local)
                    m.log_level > 99 && println("[DEBUG] VAR$(i): SOL=$(round(point,4)) RATIO=$(m.discretization_ratio)  |$(round(lb_local,2)) |$(round(lb_new,2)) <- * -> $(round(ub_new,2))| $(round(ub_local,2))|")
                    if ub_new < ub_local && !isapprox(ub_new, ub_local; atol=m.tolerance)  # Insert new UB-based partition
                        insert!(m.discretization[i], j+1, ub_new)
                    end
                    if lb_new > lb_local && !isapprox(lb_new, lb_local; atol=m.tolerance)  # Insert new LB-based partition
                        insert!(m.discretization[i], j+1, lb_new)
                    end
                    break
                end
            end
        end
    end

    return
end
