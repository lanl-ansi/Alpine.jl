"""
    init_tight_bound(m::PODNonlinearModel)

Initial tight bound vectors for later operations
"""
function init_tight_bound(m::PODNonlinearModel)

    m.l_var_tight = [m.l_var_orig, fill(-Inf, m.num_var_linear_mip+m.num_var_nonlinear_mip);]
    m.u_var_tight = [m.u_var_orig, fill(Inf, m.num_var_linear_mip+m.num_var_nonlinear_mip);]
    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Bin
            m.l_var_tight[i] = 0.0
            m.u_var_tight[i] = 1.0
        elseif m.var_type_orig[i] == :Int
            m.l_var_tight[i] = floor(m.l_var_tight[i])
            m.u_var_tight[i] = ceil(m.u_var_tight[i])
        end
    end

    return
end

"""
    init_disc(m::PODNonlinearModel)

This function initialize the dynamic discretization used for any bounding models. By default, it takes (.l_var_orig, .u_var_orig) as the base information. User is allowed to use alternative bounds for initializing the discretization dictionary.
The output is a dictionary with MathProgBase variable indices keys attached to the :PODNonlinearModel.discretization.
"""
function init_disc(m::PODNonlinearModel)

    for var in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)
        if m.var_type[var] in [:Bin, :Cont]
            lb = m.l_var_tight[var]
            ub = m.u_var_tight[var]
            m.discretization[var] = [lb, ub]
        elseif m.var_type[var] in [:Int]
            m.int_enable ? lb = floor(m.l_var_tight[var]) - 0.5 : lb = floor(m.l_var_tight[var])
            m.int_enable ? ub = ceil(m.u_var_tight[var]) + 0.5 : ub = floor(m.u_var_tight[var])
            m.discretization[var] = [lb, ub]
        else
            error("[EXCEPTION] Unexpected variable type when initializing discretization dictionary.")
        end
    end

    return
end

"""

    to_discretization(m::PODNonlinearModel, lbs::Vector{Float64}, ubs::Vector{Float64})

Utility functions to convert bounds vectors to Dictionary based structures that is more suitable for
partition operations.

"""
function to_discretization(m::PODNonlinearModel, lbs::Vector{Float64}, ubs::Vector{Float64})

    @assert length(lbs) == length(ubs)
    var_discretization = Dict()
    for var in 1:m.num_var_orig
        lb = lbs[var]
        ub = ubs[var]
        var_discretization[var] = [lb, ub]
    end

    total_var_cnt = m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip
    orig_var_cnt = m.num_var_orig

    if length(lbs) == total_var_cnt
        for var in (1+orig_var_cnt):total_var_cnt
            lb = lbs[var]
            ub = ubs[var]
            var_discretization[var] = [lb, ub]
        end
    else
        for var in (1+orig_var_cnt):total_var_cnt
            lb = -Inf
            ub = Inf
            var_discretization[var] = [lb, ub]
        end
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
        flatten_discretization[var] = [discretization[var][1],discretization[var][end]]
    end

    return flatten_discretization
end


"""
    detect_bound_from_aff(m::PODNonlinearModel)

Detect bounds from parse affine constraint. This function examines the one variable constraints such as
x >= 5, x <= 5 or x == 5 and fetch the information to m.l_var_tight and m.u_var_tight.
This function can potential grow to be smarter.
"""
function bounds_propagation(m::PODNonlinearModel)

    exhausted = false
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
                        (m.loglevel > 99) && println("[VAR$(var_idx)] LB $(m.l_var_tight[var_idx]) evaluated from constraint")
                    end
                    if eval_u_bound < m.u_var_tight[var_idx] - m.tol
                        exhausted = false
                        m.u_var_tight[var_idx] = eval_u_bound
                        (m.loglevel > 99) && println("[VAR$(var_idx)] UB $(m.u_var_tight[var_idx]) evaluated from constraints")
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
                        (m.loglevel > 99) && println("[VAR$(var_idx)] LB $(m.l_var_tight[var_idx]) evaluated from constraints")
                    end
                elseif aff[:sense] == :(>=) && var_coef < 0.0  # -a($) - by + cz >= 100, y∈[1,10], z∈[2,50], a,b,c > 0
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
                        (m.loglevel > 99) && println("[VAR$(var_idx)] UB $(m.u_var_tight[var_idx]) evaluated from constraints")
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
                        (m.loglevel > 99) && println("[VAR$(var_idx)] UB $(m.u_var_tight[var_idx]) evaluated from constraints")
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
                        (m.loglevel > 99) && println("[VAR$(var_idx)] LB $(m.l_var_tight[var_idx]) evaluated from constraints")
                    end
                end
            end
        end
        (exhausted == true && m.loglevel > 99) && println("Initial constraint-based bound evaluation exhausted...")
    end



    return
end

"""
    Categorize variable based on variable bounds
"""
function recategorize_var(m::PODNonlinearModel)

    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Int && m.l_var_orig[i] == 0.0 && m.u_var_orig[i] == 1.0
            m.var_type_orig[i] = :Bin
            m.var_type[i] = :Bin
            println("Converting VAR$(i) to binary variable")
        end
    end

    return
end

"""
    resolve_var_bounds(m::PODNonlinearModel)

Resolve the bounds of the lifted variable using the information in l_var_tight and u_var_tight. This method only takes
in known or trivial bounds information to reason lifted variable bound to avoid the cases of infinity bounds.
"""
function resolve_var_bounds(m::PODNonlinearModel)

    # Basic Bound propagation
    m.presolve_bp && bounds_propagation(m) # Fetch bounds from constraints

    # First resolve infinity bounds with assumptions
    resolve_inf_bounds(m) # Temporarily disabled

    # Added sequential bound resolving process base on DFS process, which ensures all bounds are secured.
    # Increased complexity from linear to square but a reasonable amount
    # Potentially, additional mapping can be applied to reduce the complexity
    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            m.nonconvex_terms[k][:bound_resolver](m, k)
        elseif haskey(m.linear_terms, k)
            basic_linear_bounds(m, k)
        else
            error("[RARE] Found homeless term key $(k) during bound resolution.")
        end
    end

    return
end

"""
    Critically assumed since POD relies on finite bound to work
"""
function resolve_inf_bounds(m::PODNonlinearModel)

    warnuser = false
    infcount = 0

    # Only specify necessary bounds
    for i in m.candidate_disc_vars
        if m.l_var_tight[i] == -Inf
            warnuser = true
            m.l_var_tight[i] = -10e3
            infcount += 1
        end
        if m.u_var_tight[i] == Inf
            warnuser = true
            m.u_var_tight[i] = 10e3
            infcount +=1
        end
    end

    warnuser && println("Inf bound detected on $(infcount) variables. Initialize with value -10e4/10e4. This may affect global optimality and performance.")

    return
end

"""
    resolve_var_bounds(nonconvex_terms::Dict, discretization::Dict)

    For discretization to be performed, we do not allow for a variable being discretized to have infinite bounds.
    The lifted variables will have infinite bounds and the function infers bounds on these variables. This process
    can help speed up the subsequent solve in subsequent iterations.

    Only used in presolve bound tightening
"""
function resolve_var_bounds(m::PODNonlinearModel, d::Dict; kwargs...)

    # Added sequential bound resolving process base on DFS process, which ensures all bounds are secured.
    # Increased complexity from linear to square but a reasonable amount
    # Potentially, additional mapping can be applied to reduce the complexity
    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            nlk = k
            if m.nonconvex_terms[nlk][:nonlinear_type] in POD_C_MONOMIAL
                d = basic_monomial_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:BININT]
                d = basic_binint_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:BINPROD]
                d = basic_binprod_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in POD_C_TRIGONOMETRIC
                d = basic_sincos_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:BINLIN]
                d = basic_binlin_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:INTLIN]
                d = basic_intlin_bounds(m, nlk, d)
            elseif m.nonconvex_terms[nlk][:nonlinear_type] in [:INTPROD]
                d = basic_intprod_bounds(m, nlk, d)
            else
                error("EXPECTED ERROR : NEED IMPLEMENTATION")
            end
        elseif haskey(m.linear_terms, k)
            d = basic_linear_bounds(m, k, d)
        else
            error("[RARE] Found homeless term key $(k) during bound resolution.")
        end
    end

    return d
end

"""
    resolve_closed_var_bounds(m::PODNonlinearModel)

This function seeks variable with tight bounds (by presolve_bt_width_tol) by checking .l_var_tight and .u_var_tight.
If a variable is found to be within a sufficiently small interval then no discretization will be performed on this variable
and the .discretization will be cleared with the tight bounds for basic McCormick operation if necessary.

"""
function resolve_closed_var_bounds(m::PODNonlinearModel; kwargs...)

    for var in m.candidate_disc_vars
        if abs(m.l_var_tight[var] - m.u_var_tight[var]) < m.presolve_bt_width_tol         # Closed Bound Criteria
            deleteat!(m.disc_vars, findfirst(m.disc_vars, var)) # Clean nonconvex_terms by deleting the info
            m.discretization[var] = [m.l_var_tight[var], m.u_var_tight[var]]              # Clean up the discretization for basic McCormick if necessary
        end
    end

    return
end

"""
    update_var_bounds(m::PODNonlinearModel, discretization::Dict; len::Float64=length(keys(discretization)))

This function take in a dictionary-based discretization information and convert them into two bounds vectors (l_var, u_var) by picking the smallest and largest numbers. User can specify a certain length that may contains variables that is out of the scope of discretization.

Output::

    l_var::Vector{Float64}, u_var::Vector{Float64}
"""
function update_var_bounds(discretization; kwargs...)

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
