"""
    initialize_tight_bounds(m::PODNonlinearModel)

Initial tight bound vectors for later operations
"""
function initialize_tight_bounds(m::PODNonlinearModel)

    m.l_var_tight = [m.l_var_orig,fill(-Inf, m.num_var_lifted_mip);]
    m.u_var_tight = [m.u_var_orig,fill(Inf, m.num_var_lifted_mip);]
    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Bin
            m.l_var_tight[i] = 0.0
            m.u_var_tight[i] = 1.0
        end
    end

    return
end


"""
    detect_bound_from_aff(m::PODNonlinearModel)

Detect bounds from parse affine constraint. This function examines the one variable constraints such as
x >= 5, x <= 5 or x == 5 and fetch the information to m.l_var_tight and m.u_var_tight.
This function can potential grow to be smarter.
"""
function detect_bound_from_aff(m::PODNonlinearModel)

    for aff in m.lifted_constr_aff_mip
        if length(aff[:coefs]) == 1 && length(aff[:vars]) == 1
            var_idx = aff[:vars][1].args[2]
            @assert (isa(var_idx, Float64) || isa(var_idx, Int))
            if aff[:sense] == :(==)
                m.l_var_tight[var_idx] = aff[:coefs][1] * aff[:rhs]
                m.u_var_tight[var_idx] = aff[:coefs][1] * aff[:rhs]
                (m.log_level > 99) && println("[VAR$(var_idx)] Both bound $(aff[:rhs]) detected from constraints")
            elseif aff[:sense] == :(>=) && aff[:coefs][1] > 0.0
                eval_bound = aff[:rhs]/aff[:coefs][1]
                if eval_bound > m.l_var_tight[var_idx] + m.tol
                    m.l_var_tight[var_idx] = eval_bound
                    (m.log_level > 99) && println("[VAR$(var_idx)] Lower bound $(m.l_var_tight[var_idx]) detected from constraints")
                end
            elseif aff[:sense] == :(>=) && aff[:coefs][1] < 0.0  # Flip sign
                eval_bound = aff[:rhs]/aff[:coefs][1]
                if eval_bound < m.u_var_tight[var_idx] - m.tol
                    m.u_var_tight[var_idx] = eval_bound
                    (m.log_level > 99) && println("[VAR$(var_idx)] Upper bound $(m.u_var_tight[var_idx]) detected from constraints")
                end
            elseif aff[:sense] == :(<=) && aff[:coefs][1] > 0.0
                eval_bound = aff[:rhs]/aff[:coefs][1]
                if eval_bound < m.u_var_tight[var_idx] - m.tol
                    m.u_var_tight[var_idx] = eval_bound
                    (m.log_level > 99) && println("[VAR$(var_idx)] Upper bound $(m.u_var_tight[var_idx]) detected from constraints")
                end
            elseif aff[:sense] == :(<=) && aff[:coefs][1] < 0.0  # Flip sign
                eval_bound = aff[:rhs]/aff[:coefs][1]
                if eval_bound > m.l_var_tight[var_idx] + m.tol
                    m.l_var_tight[var_idx] = eval_bound
                    (m.log_level > 99) && println("[VAR$(var_idx)] Lower bound $(m.l_var_tight[var_idx]) detected from constraints")
                end
            end
        end
    end

    exhausted = false
    while !exhausted
        exhausted = true
        for aff in m.lifted_constr_aff_mip
            for i in 1:length(aff[:vars])
                var_idx = aff[:vars][i].args[2]
                var_coef = aff[:coefs][i]
                @assert (isa(var_idx, Float64) || isa(var_idx, Int))
                if aff[:sense] == :(==) && var_coef > 0.0
                    eval_l_bound = aff[:rhs] / var_coef
                    eval_u_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j]*var_coef > 0.0     # same sign
                            (eval_l_bound != -Inf) && (eval_l_bound -= abs(aff[:coefs][j]/var_coef)*m.u_var_tight[aff[:vars][j].args[2]])
                            (eval_u_bound != Inf) && (eval_u_bound -= abs(aff[:coefs][j]/var_coef)*m.l_var_tight[aff[:vars][j].args[2]])
                        elseif j!= i && aff[:coefs][j]*var_coef < 0.0  # different sign
                            (eval_l_bound != -Inf) && (eval_l_bound += abs(aff[:coefs][j]/var_coef)*m.l_var_tight[aff[:vars][j].args[2]])
                            (eval_u_bound != Inf) && (eval_u_bound += abs(aff[:coefs][j]/var_coef)*m.u_var_tight[aff[:vars][j].args[2]])
                        end
                    end
                    if eval_l_bound > m.l_var_tight[var_idx] + m.tol
                        exhausted = false
                        m.l_var_tight[var_idx] = eval_l_bound
                        (m.log_level > 99) && println("[VAR$(var_idx)] Lower bound $(m.l_var_tight[var_idx]) evaluated from constraints")
                    end
                    if eval_u_bound < m.u_var_tight[var_idx] - m.tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_u_bound
                        (m.log_level > 99) && println("[VAR$(var_idx)] Upper bound $(m.u_var_tight[var_idx]) evaluated from constraints")
                    end
                elseif aff[:sense] == :(>=) && var_coef > 0.0  # a($) - by + cz >= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound -= abs(aff[:coefs][j]/var_coef) * m.u_var_tight[aff[:vars][j].args[2]]
                        elseif j !=i aff[:coefs][j] < 0.0
                            eval_bound += abs(aff[:coefs][j]/var_coef) * m.l_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == -Inf) && break
                    end
                    if eval_bound > m.l_var_tight[var_idx] + m.tol
                        exhausted = false
                        m.l_var_tight[var_idx] = eval_bound
                        (m.log_level > 99) && println("[VAR$(var_idx)] Lower bound $(m.l_var_tight[var_idx]) evaluated from constraints")
                    end
                elseif var_coef < 0.0 && aff[:sense] == :(>=) # -a($) - by + cz >= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound += abs(aff[:coefs][j]/var_coef) * m.u_var_tight[aff[:vars][j].args[2]]
                        elseif j != i && aff[:coefs][j] < 0.0
                            eval_bound -= abs(aff[:coefs][j]/var_coef) * m.l_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == Inf) && break
                    end
                    if eval_bound < m.u_var_tight[var_idx] - m.tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_bound
                        (m.log_level > 99) && println("[VAR$(var_idx)] Upper bound $(m.u_var_tight[var_idx]) evaluated from constraints")
                    end
                elseif (aff[:sense] == :(<=) && aff[:coefs][i] > 0.0) # a($) - by + cz <= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound -= abs(aff[:coefs][j]/var_coef) * m.l_var_tight[aff[:vars][j].args[2]]
                        elseif j != i && aff[:coefs][j] < 0.0
                            eval_bound += abs(aff[:coefs][j]/var_coef) * m.u_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == Inf) && break
                    end
                    if eval_bound < m.u_var_tight[var_idx] - m.tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_bound
                        (m.log_level > 99) && println("[VAR$(var_idx)] Upper bound $(m.u_var_tight[var_idx]) evaluated from constraints")
                    end
                elseif (aff[:sense] == :(<=) && aff[:coefs][i] < 0.0) # -a($) - by + cz <= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound += abs(aff[:coefs][j]/var_coef) * m.l_var_tight[aff[:vars][j].args[2]]
                        elseif j != i && aff[:coefs][j] < 0.0
                            eval_bound -= abs(aff[:coefs][j]/var_coef) * m.u_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == -Inf) && break
                    end
                    if eval_bound > m.l_var_tight[var_idx] + m.tol
                        exhausted = false
                        m.l_var_tight[var_idx] = eval_bound
                        (m.log_level > 99) && println("[VAR$(var_idx)] Lower bound $(m.l_var_tight[var_idx]) evaluated from constraints")
                    end
                end
            end
        end
        (exhausted == true && m.log_level > 99) && println("Initial constraint-based bound evaluation exhausted...")
    end

    return
end


"""
    resolve_lifted_var_bounds(m::PODNonlinearModel)

Resolve the bounds of the lifted variable using the information in l_var_tight and u_var_tight. This method only takes
in known or trivial bounds information to reason lifted variable bound to avoid the cases of infinity bounds.
"""
function resolve_lifted_var_bounds(m::PODNonlinearModel)

    # Added sequential bound resolving process base on DFS process, which ensures all bounds are secured.
    # Increased complexity from linear to square but a reasonable amount
    # Potentially, additional mapping can be applied to reduce the complexity
    for i in 1:length(m.nonlinear_info)
        for bi in keys(m.nonlinear_info)
            if (m.nonlinear_info[bi][:id]) == i && (m.nonlinear_info[bi][:nonlinear_type] in [:bilinear, :monomial])
                idx_a = bi[1].args[2]
                idx_b = bi[2].args[2]
                idx_ab = m.nonlinear_info[bi][:lifted_var_ref].args[2]
                bound = [m.l_var_tight[idx_a], m.u_var_tight[idx_a]] * [m.l_var_tight[idx_b], m.u_var_tight[idx_b]]'
                if minimum(bound) > m.l_var_tight[idx_ab] + m.tol
                    m.l_var_tight[idx_ab] = minimum(bound)
                end
                if maximum(bound) < m.u_var_tight[idx_ab] + m.tol
                    m.u_var_tight[idx_ab] = maximum(bound)
                end
            end
        end
    end

    return
end

"""
    resolve_lifted_var_bounds(nonlinear_info::Dict, discretization::Dict)

For discretization to be performed, we do not allow for a variable being discretized to have infinite bounds.
The lifted variables will have infinite bounds and the function infers bounds on these variables. This process
can help speed up the subsequent solve in subsequent iterations.
"""
function resolve_lifted_var_bounds(nonlinear_info::Dict, discretization::Dict; kwargs...)

    # Added sequential bound resolving process base on DFS process, which ensures all bounds are secured.
    # Increased complexity from linear to square but a reasonable amount
    # Potentially, additional mapping can be applied to reduce the complexity
    for i in 1:length(nonlinear_info)
        for bi in keys(nonlinear_info)
            if (nonlinear_info[bi][:id] == i) && (nonlinear_info[bi][:nonlinear_type] in [:bilinear, :monomial])
                idx_a = bi[1].args[2]
                idx_b = bi[2].args[2]
                idx_ab = nonlinear_info[bi][:lifted_var_ref].args[2]
                bound = [discretization[idx_a][1], discretization[idx_a][end]] * [discretization[idx_b][1], discretization[idx_b][end]]'
                discretization[idx_ab] = [-Inf, Inf]
                discretization[idx_ab][1] = minimum(bound)
                discretization[idx_ab][2] = maximum(bound)
            end
        end
    end

    return discretization
end

"""
    resolve_closed_var_bounds(m::PODNonlinearModel)

This function seeks variable with tight bounds (by presolve_bt_width_tol) by checking .l_var_tight and .u_var_tight.
If a variable is found to be within a sufficiently small interval then no discretization will be performed on this variable
and the .discretization will be cleared with the tight bounds for basic McCormick operation if necessary.

"""
function resolve_closed_var_bounds(m::PODNonlinearModel; kwargs...)

    for var in m.all_nonlinear_vars
        if abs(m.l_var_tight[var] - m.u_var_tight[var]) < m.presolve_bt_width_tol         # Closed Bound Criteria
            deleteat!(m.var_discretization_mip, findfirst(m.var_discretization_mip, var)) # Clean nonlinear_info by deleting the info
            m.discretization[var] = [m.l_var_tight[var], m.u_var_tight[var]]              # Clean up the discretization for basic McCormick if necessary
        end
    end

    return
end
