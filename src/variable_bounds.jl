"""
    init_tight_bound(m::Optimizer)

Initialize internal bound vectors (placeholders) to be used in other places.
In this case, we don't have to mess with the original bound information.
"""
function init_tight_bound(m::Optimizer)
    m.l_var_tight = [
        m.l_var_orig
        fill(-Inf, m.num_var_linear_mip + m.num_var_nonlinear_mip)
    ]
    m.u_var_tight = [
        m.u_var_orig
        fill(Inf, m.num_var_linear_mip + m.num_var_nonlinear_mip)
    ]
    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Bin
            m.l_var_tight[i] = 0.0
            m.u_var_tight[i] = 1.0
            # elseif m.var_type_orig[i] == :Int
            #     m.l_var_tight[i] = floor(m.l_var_tight[i])
            #     m.u_var_tight[i] = ceil(m.u_var_tight[i])
        end
    end

    return
end

"""
    init_disc(m::Optimizer)

This function initialize the dynamic discretization used for any bounding models. By default, it takes (.l_var_orig, .u_var_orig) as the base information. User is allowed to use alternative bounds for initializing the discretization dictionary.
The output is a dictionary with MathProgBase variable indices keys attached to the :Optimizer.discretization.
"""
function init_disc(m::Optimizer)
    for var in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)
        if m.var_type[var] in [:Bin, :Cont]
            lb = m.l_var_tight[var]
            ub = m.u_var_tight[var]
            m.discretization[var] = [lb, ub]
        else
            error(
                "[EXCEPTION] Unexpected variable type when initializing discretization dictionary.",
            )
        end
    end

    return
end

"""

    _get_discretization_dict(m::Optimizer, lbs::Vector{Float64}, ubs::Vector{Float64})

Utility functions to convert bounds vectors to Dictionary based structures that are more suitable for
partition operations.

"""
function _get_discretization_dict(
    m::Optimizer,
    lbs::Vector{Float64},
    ubs::Vector{Float64},
)
    @assert length(lbs) == length(ubs)

    var_discretization = Dict()
    total_var_cnt = m.num_var_orig + m.num_var_linear_mip + m.num_var_nonlinear_mip

    for var in 1:total_var_cnt
        if length(lbs) == total_var_cnt
            lb = lbs[var]
            ub = ubs[var]
        else
            lb = -Inf
            ub = Inf
        end
        var_discretization[var] = [lb, ub]
    end

    return var_discretization
end

"""
    flatten_discretization(discretization::Dict)

Utility functions to eliminate all partition on discretizing variable and keep the loose bounds.

"""
function flatten_discretization(discretization::Dict; kwargs...)
    flatten_discretization = Dict()
    for var in keys(discretization)
        flatten_discretization[var] = [discretization[var][1], discretization[var][end]]
    end

    return flatten_discretization
end

"""
    detect_bound_from_aff(m::Optimizer)

Detect bounds from parse affine constraint. This function examines the one variable constraints such as
x >= 5, x <= 5 or x == 5 and fetch the information to m.l_var_tight and m.u_var_tight.
"""
function bound_propagation(m::Optimizer)
    exhausted = false
    infeasible = false
    tol = Alp.get_option(m, :tol)

    while !exhausted
        exhausted = true
        for aff in m.bounding_constr_mip
            for i in 1:length(aff[:vars])
                var_idx = aff[:vars][i].args[2]
                var_coef = aff[:coefs][i]
                @assert (isa(var_idx, Float64) || isa(var_idx, Int))

                if aff[:sense] == :(==) && var_coef > 0.0
                    eval_l_bound = aff[:rhs] / var_coef
                    eval_u_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] * var_coef > 0.0     # same sign
                            (eval_l_bound != -Inf) && (
                                eval_l_bound -=
                                    abs(aff[:coefs][j] / var_coef) *
                                    m.u_var_tight[aff[:vars][j].args[2]]
                            )
                            (eval_u_bound != Inf) && (
                                eval_u_bound -=
                                    abs(aff[:coefs][j] / var_coef) *
                                    m.l_var_tight[aff[:vars][j].args[2]]
                            )
                        elseif j != i && aff[:coefs][j] * var_coef < 0.0  # different sign
                            (eval_l_bound != -Inf) && (
                                eval_l_bound +=
                                    abs(aff[:coefs][j] / var_coef) *
                                    m.l_var_tight[aff[:vars][j].args[2]]
                            )
                            (eval_u_bound != Inf) && (
                                eval_u_bound +=
                                    abs(aff[:coefs][j] / var_coef) *
                                    m.u_var_tight[aff[:vars][j].args[2]]
                            )
                        end
                    end

                    if eval_l_bound > m.l_var_tight[var_idx] + tol
                        exhausted = false
                        m.l_var_tight[var_idx] = eval_l_bound
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] LB $(m.l_var_tight[var_idx]) evaluated from constraint",
                        )
                    elseif eval_l_bound > m.u_var_tight[var_idx] + tol
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] Infeasibility detection during bound propagation",
                        )
                        infeasible = true
                        break
                    end

                    if eval_u_bound < m.u_var_tight[var_idx] - tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_u_bound
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] UB $(m.u_var_tight[var_idx]) evaluated from constraints",
                        )
                    elseif eval_u_bound < m.l_var_tight[var_idx] - tol
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] Infeasibility detection during bound propagation",
                        )
                        infeasible = true
                        break
                    end

                elseif aff[:sense] == :(>=) && var_coef > 0.0  # a($) - by + cz >= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound -=
                                abs(aff[:coefs][j] / var_coef) *
                                m.u_var_tight[aff[:vars][j].args[2]]
                        elseif j != i
                            aff[:coefs][j] < 0.0
                            eval_bound +=
                                abs(aff[:coefs][j] / var_coef) *
                                m.l_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == -Inf) && break
                    end

                    if eval_bound > m.l_var_tight[var_idx] + tol
                        exhausted = false
                        m.l_var_tight[var_idx] = eval_bound
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] LB $(m.l_var_tight[var_idx]) evaluated from constraints",
                        )
                    elseif eval_bound > m.u_var_tight[var_idx] + tol
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] Infeasibility detection during bound propagation",
                        )
                        infeasible = true
                        break
                    end

                elseif aff[:sense] == :(>=) && var_coef < 0.0  # -a($) - by + cz >= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound +=
                                abs(aff[:coefs][j] / var_coef) *
                                m.u_var_tight[aff[:vars][j].args[2]]
                        elseif j != i && aff[:coefs][j] < 0.0
                            eval_bound -=
                                abs(aff[:coefs][j] / var_coef) *
                                m.l_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == Inf) && break
                    end

                    if eval_bound < m.u_var_tight[var_idx] - tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_bound
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] UB $(m.u_var_tight[var_idx]) evaluated from constraints",
                        )
                    elseif eval_bound < m.l_var_tight[var_idx] - tol
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] Infeasibility detection during bound propagation",
                        )
                        infeasible = true
                        break
                    end

                elseif (aff[:sense] == :(<=) && aff[:coefs][i] > 0.0) # a($) - by + cz <= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound -=
                                abs(aff[:coefs][j] / var_coef) *
                                m.l_var_tight[aff[:vars][j].args[2]]
                        elseif j != i && aff[:coefs][j] < 0.0
                            eval_bound +=
                                abs(aff[:coefs][j] / var_coef) *
                                m.u_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == Inf) && break
                    end

                    if eval_bound < m.u_var_tight[var_idx] - tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_bound
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] UB $(m.u_var_tight[var_idx]) evaluated from constraints",
                        )
                    elseif eval_bound < m.l_var_tight[var_idx] - tol
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] Infeasibility detection during bound propagation",
                        )
                        infeasible = true
                        break
                    end
                elseif (aff[:sense] == :(<=) && aff[:coefs][i] < 0.0) # -a($) - by + cz <= 100, y∈[1,10], z∈[2,50], a,b,c > 0
                    eval_bound = aff[:rhs] / var_coef
                    for j in 1:length(aff[:vars])
                        if j != i && aff[:coefs][j] > 0.0
                            eval_bound +=
                                abs(aff[:coefs][j] / var_coef) *
                                m.l_var_tight[aff[:vars][j].args[2]]
                        elseif j != i && aff[:coefs][j] < 0.0
                            eval_bound -=
                                abs(aff[:coefs][j] / var_coef) *
                                m.u_var_tight[aff[:vars][j].args[2]]
                        end
                        (eval_bound == -Inf) && break
                    end

                    if eval_bound > m.l_var_tight[var_idx] + tol
                        exhausted = false
                        m.l_var_tight[var_idx] = eval_bound
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] LB $(m.l_var_tight[var_idx]) evaluated from constraints",
                        )
                    elseif eval_bound > m.u_var_tight[var_idx] + tol
                        (Alp.get_option(m, :log_level) > 199) && println(
                            "[VAR$(var_idx)] Infeasibility detection during bound propagation",
                        )
                        infeasible = true
                        break
                    end
                end
            end
        end
        (exhausted == true && Alp.get_option(m, :log_level) > 99) &&
            println("Initial constraint-based bound evaluation exhausted...")
    end

    if infeasible
        m.status[:bounding_solve] = MOI.INFEASIBLE
        @warn "[INFEASIBLE] Infeasibility detected via bound propagation"
    end

    return infeasible
end

"""
    Recategorize :Int variables to :Bin variables if variable bounds are [0,1]
"""
function recategorize_var(m::Optimizer)
    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Int && m.l_var_orig[i] == 0.0 && m.u_var_orig[i] == 1.0
            m.var_type_orig[i] = :Bin
            m.var_type[i] = :Bin
        end
    end

    return
end

"""
    resolve_var_bounds(m::Optimizer)

Resolve the bounds of the lifted variable using the information in l_var_tight and u_var_tight. This method only takes
in known or trivial bounds information to deduce lifted variable bounds and to potentially avoid the cases of infinity bounds.
"""
function resolve_var_bounds(m::Optimizer)

    # Basic Bound propagation
    if Alp.get_option(m, :presolve_bp)
        setproperty!(m, :presolve_infeasible, Alp.bound_propagation(m)) # Fetch bounds from constraints
    end

    # Resolve unbounded variables in the original formulation 
    Alp.resolve_inf_bounds(m)

    # Added sequential bound resolving process base on DFS process, which ensures all bounds are secured.
    # Increased complexity from linear to square but a reasonable amount
    # Potentially, additional mapping can be applied to reduce the complexity
    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            m.nonconvex_terms[k][:bound_resolver](m, k)
        elseif haskey(m.linear_terms, k)
            Alp.basic_linear_bounds(m, k)
        else
            error("Found homeless term key $(k) during bound resolution.")
        end
    end

    return
end

"""
    Resolving Inf variable bounds since Alpine relies on finite bounds to construct relaxations
"""
function resolve_inf_bounds(m::Optimizer)
    warnuser = false
    infcount_l = 0
    infcount_u = 0

    # Only specify necessary bounds
    for i in 1:length(m.l_var_orig)
        if m.l_var_tight[i] == -Inf
            warnuser = true
            m.l_var_tight[i] = -Alp.get_option(m, :large_bound)
            infcount_l += 1
        end
        if m.u_var_tight[i] == Inf
            warnuser = true
            m.u_var_tight[i] = Alp.get_option(m, :large_bound)
            infcount_u += 1
        end
    end
    infcount = min(infcount_l, infcount_u)
    if infcount == 1
        warnuser && println(
            "Warning: -/+Inf bounds detected on at least $infcount variable. Initializing with values -/+$(Alp.get_option(m, :large_bound)). This may affect global optimal values and run times.",
        )
    elseif infcount > 1
        warnuser && println(
            "Warning: -/+Inf bounds detected on at least $infcount variables. Initializing with values -/+$(Alp.get_option(m, :large_bound)). This may affect global optimal values and run times.",
        )
    end

    return
end

"""
    resolve_var_bounds(nonconvex_terms::Dict, discretization::Dict)

    For discretization to be performed, we do not allow a variable being discretized to have infinite bounds.
    The lifted/auxiliary variables may have infinite bounds and the function infers bounds on these variables. This process
    can help speed up the subsequent solve times.

    Only used in presolve bound tightening
"""
function resolve_var_bounds(m::Optimizer, d::Dict; kwargs...)
    # Added sequential bound resolving process base on DFS process, which ensures all bounds are accurate.
    # Potentially, additional mapping can be applied to reduce the complexity
    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            nlk = k
            if m.nonconvex_terms[nlk][:nonlinear_type] in ALPINE_C_MONOMIAL
                d = Alp.basic_monomial_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:BINPROD]
                d = Alp.basic_binprod_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:BINLIN]
                d = Alp.basic_binlin_bounds(m, nlk, d)
            else
                error("EXPECTED ERROR : NEED IMPLEMENTATION")
            end
        elseif haskey(m.linear_terms, k)
            d = Alp.basic_linear_bounds(m, k, d)
        else
            error(
                "Found a homeless term key $(k) during bound resolution im resolve_var_bounds.",
            )
        end
    end

    return d
end

"""
    update_var_bounds(m::Optimizer, discretization::Dict; len::Float64=length(keys(discretization)))

This function takes in a dictionary-based discretization information and convert them into two bounds vectors (l_var, u_var) by picking the smallest and largest numbers. User can specify a certain length that may contains variables that is out of the scope of discretization.

Output::

    l_var::Vector{Float64}, u_var::Vector{Float64}
"""
function update_var_bounds(discretization::Dict{Any,Any}; kwargs...)
    options = Dict(kwargs)

    haskey(options, :len) ? len = options[:len] : len = length(keys(discretization))

    l_var = fill(-Inf, len)
    u_var = fill(Inf, len)

    for var_idx in keys(discretization)
        l_var[var_idx] = discretization[var_idx][1]
        u_var[var_idx] = discretization[var_idx][end]
    end

    return l_var, u_var
end

# ================================================================= unused, but needs testing
# """
#     resolve_closed_var_bounds(m::Optimizer)

# This function seeks variable with tight bounds (by presolve_bt_width_tol) by checking .l_var_tight and .u_var_tight.
# If a variable is found to be within a sufficiently small interval then no discretization will be performed on this variable
# and the .discretization will be cleared with the tight bounds for basic McCormick operation if necessary.

# """
# function resolve_closed_var_bounds(m::Optimizer; kwargs...)
#     for var in m.candidate_disc_vars
#         if abs(m.l_var_tight[var] - m.u_var_tight[var]) <
#            Alp.get_option(m, :presolve_bt_width_tol)         # Closed Bound Criteria
#             deleteat!(m.disc_vars, findfirst(m.disc_vars, var)) # Clean nonconvex_terms by deleting the info
#             m.discretization[var] = [m.l_var_tight[var], m.u_var_tight[var]]              # Clean up the discretization for basic McCormick if necessary
#         end
#     end

#     return
# end
