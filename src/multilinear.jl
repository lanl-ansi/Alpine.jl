function amp_post_convhull(m::Optimizer; kwargs...)
    options = Dict(kwargs)
    haskey(options, :use_disc) ? d = options[:use_disc] : d = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables
    β = Dict()  # Lifted variables for exact formulation

    # Convexification for non-convex terms
    contains_multilinear = false
    for k in keys(m.nonconvex_terms)
        nl_type = m.nonconvex_terms[k][:nonlinear_type]
        if ((nl_type == :MULTILINEAR) || (nl_type == :BILINEAR)) &&
           (m.nonconvex_terms[k][:convexified] == false)
            λ, α = Alp.amp_convexify_multilinear(m, k, λ, α, d)
            contains_multilinear = true
        elseif nl_type == :MONOMIAL && !m.nonconvex_terms[k][:convexified]
            λ, α = Alp.amp_convexify_quadratic_univariate(m, k, λ, α, d)
        elseif nl_type == :BINLIN && !m.nonconvex_terms[k][:convexified]
            β = Alp.amp_convexify_binlin(m, k, β)
        elseif nl_type == :BINPROD && !m.nonconvex_terms[k][:convexified]
            β = Alp.amp_convexify_multilinear_binary(m, k, β)
        end
    end

    # Add lambda linking constraints
    if m.options.linking_constraints && contains_multilinear
        Alp._add_multilinear_linking_constraints(m, λ)
    end

    # Code for warm starting bounding MIP iterations
    Alp.get_option(m, :convhull_warmstart) &&
        !isempty(m.best_bound_sol) &&
        Alp.amp_warmstart_α(m, α)

    return
end

function amp_convexify_multilinear(
    m::Optimizer,
    k::Any,
    λ::Dict,
    α::Dict,
    discretization::Dict,
)
    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    ml_indices, dim, extreme_point_cnt = Alp.amp_convhull_prepare(m, discretization, k)   # convert key to easy read mode
    λ = Alp.amp_convhull_λ(m, k, ml_indices, λ, extreme_point_cnt, dim)
    λ = Alp.populate_convhull_extreme_values(
        m,
        discretization,
        ml_indices,
        λ,
        dim,
        ones(Int, length(dim)),
    )
    α = Alp.amp_convhull_α(m, ml_indices, α, dim, discretization)
    Alp.amp_post_convhull_constrs(
        m,
        λ,
        α,
        ml_indices,
        dim,
        extreme_point_cnt,
        discretization,
    )

    return λ, α
end

function amp_convexify_quadratic_univariate(
    m::Optimizer,
    k::Any,
    λ::Dict,
    α::Dict,
    discretization::Dict,
)
    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    monomial_index, dim, extreme_point_cnt =
        Alp.amp_convhull_prepare(m, discretization, k, monomial = true)
    λ = Alp.amp_convhull_λ(m, k, monomial_index, λ, extreme_point_cnt, dim)
    λ = Alp.populate_convhull_extreme_values(m, discretization, monomial_index, λ, 2)
    α = Alp.amp_convhull_α(m, [monomial_index], α, dim, discretization)
    Alp.amp_post_convhull_constrs(m, λ, α, monomial_index, discretization)

    return λ, α
end

function amp_convexify_binlin(m::Optimizer, k::Any, β::Dict)
    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    @assert length(m.nonconvex_terms[k][:var_idxs]) == 2

    lift_idx = m.nonconvex_terms[k][:y_idx]

    if haskey(β, lift_idx)
        return β
    else
        β[lift_idx] = _index_to_variable_ref(m.model_mip, lift_idx)
    end

    bin_idx = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Bin]
    cont_idx = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] != :Bin]

    @assert length(bin_idx) == length(cont_idx) == 1

    bin_idx = bin_idx[1]
    cont_idx = cont_idx[1]

    Alp.relaxation_bilinear(
        m.model_mip,
        _index_to_variable_ref(m.model_mip, lift_idx),
        _index_to_variable_ref(m.model_mip, bin_idx),
        _index_to_variable_ref(m.model_mip, cont_idx),
        0,
        1,
        m.l_var_tight[cont_idx],
        m.u_var_tight[cont_idx],
    )

    return β
end

function amp_convexify_multilinear_binary(m::Optimizer, k::Any, β::Dict)
    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    lift_idx = m.nonconvex_terms[k][:y_idx]
    if haskey(β, lift_idx)
        return β    # Already constructed
    else
        β[lift_idx] = _index_to_variable_ref(m.model_mip, lift_idx)
    end

    z = _index_to_variable_ref(m.model_mip, m.nonconvex_terms[k][:y_idx])
    x = [_index_to_variable_ref(m.model_mip, i) for i in m.nonconvex_terms[k][:var_idxs]]

    if length(x) >= 1
        Alp.relaxation_multilinear_binary(m.model_mip, z, x)
    end

    return β
end

"""
    Method for general nonlinear terms
"""
function amp_convhull_prepare(m::Optimizer, d::Dict, nonlinear_key::Any; monomial = false)
    counted_var = []                # Keep both vector and set for collection sake
    id = Set()                      # Coverting the nonlinear indices into a set

    if isa(nonlinear_key, Vector)
        for var in nonlinear_key        # This output regulates the sequence of how composing variable should be arranged
            @assert isa(var, Expr)
            m.var_type[var.args[2]] in [:Cont, :Int] && push!(id, var.args[2])
            m.var_type[var.args[2]] in [:Cont, :Int] && push!(counted_var, var.args[2])
        end
    elseif isa(nonlinear_key, Dict)
        for var in nonlinear_key[:vars]
            @assert isa(var, Int)
            m.var_type[var] in [:Cont, :Int] && push!(id, var)
            m.var_type[var] in [:Cont, :Int] && push!(counted_var, var)
        end
    end

    if length(id) < length(counted_var) # Got repeating terms, now the sequence matters
        id = []
        for var in nonlinear_key
            m.var_type[var.args[2]] in [:Cont, :Int] && push!(id, var.args[2])
        end
    end

    dim = [length(d[i]) for i in id]

    monomial && return id[1], tuple(dim[1]), dim[1]   # One less dimension is required
    return id, tuple([i for i in dim]...), prod(dim)
end

"""
    Method for integers
"""
function amp_convhull_prepare(m::Optimizer, d::Dict, idx::Int)
    return [idx], tuple(length(d[idx])), length(d[idx])
end

"""
    Method for general nonlinear terms
"""
function amp_convhull_λ(
    m::Optimizer,
    nonlinear_key::Any,
    indices::Any,
    λ::Dict,
    ext_cnt::Int,
    dim::Tuple,
)
    y_idx = m.nonconvex_terms[nonlinear_key][:y_idx]

    @assert !(y_idx in keys(λ))
    λ[indices] = Dict(
        :dim => dim,
        :lifted_var_idx => y_idx,
        :indices => reshape([1:ext_cnt;], dim),
        :vars => JuMP.@variable(
            m.model_mip,
            [1:ext_cnt],
            lower_bound = 0,
            base_name = "L$(y_idx)"
        ),
        :vals => ones(dim),
    )

    return λ
end

"""
    Method for power terms
"""
function populate_convhull_extreme_values(
    m::Optimizer,
    d::Dict,
    mono_idx::Int,
    λ::Dict,
    p::Int,
)
    λ[mono_idx][:vals] = [d[mono_idx][i]^p for i in 1:length(d[mono_idx])]
    return λ
end

"""
    Method for regular multilinear terms
"""
function populate_convhull_extreme_values(
    m::Optimizer,
    discretization::Dict,
    indices::Any,
    λ::Dict,
    dim::Tuple,
    locator::Array,
    level::Int = 1,
)
    if level > length(dim)
        @assert length(indices) == length(dim)
        @assert length(indices) == length(locator)
        val = 1.0
        k = 0
        for i in indices
            k += 1
            val *= discretization[i][locator[k]]                      # Calculate extreme point z-value
        end
        λ[indices][:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i
            λ = Alp.populate_convhull_extreme_values(
                m,
                discretization,
                indices,
                λ,
                dim,
                locator,
                level + 1,
            )
        end
    end

    return λ
end

"""
    General Method for all term
"""
function amp_convhull_α(
    m::Optimizer,
    indices::Any,
    α::Dict,
    dim::Tuple,
    discretization::Dict,
)
    for i in indices
        if !(i in keys(α))
            lambda_cnt = length(discretization[i])
            partition_cnt = length(discretization[i]) - 1
            if Alp.get_option(m, :convhull_ebd) && partition_cnt > 2
                αCnt = Int(ceil(log(2, partition_cnt)))
                α[i] = JuMP.@variable(
                    m.model_mip,
                    [1:αCnt],
                    Bin,
                    base_name = string("YL", i)
                )
            else
                α[i] = JuMP.@variable(
                    m.model_mip,
                    [1:partition_cnt],
                    Bin,
                    base_name = "A$(i)"
                )
                JuMP.@constraint(m.model_mip, sum(α[i]) == 1)
                JuMP.@constraint(
                    m.model_mip,
                    _index_to_variable_ref(m.model_mip, i) >=
                    sum(α[i][j] * discretization[i][j] for j in 1:lambda_cnt-1)
                ) # Add x = f(α) for regulating the domains
                JuMP.@constraint(
                    m.model_mip,
                    _index_to_variable_ref(m.model_mip, i) <=
                    sum(α[i][j-1] * discretization[i][j] for j in 2:lambda_cnt)
                )
            end
        end
    end

    return α
end

function amp_convhull_α(m::Optimizer, idx::Int, α::Dict, dim, d::Dict)
    return amp_convhull_α(m, [idx], α, dim, d)
end

function amp_warmstart_α(m::Optimizer, α::Dict)
    d = m.discretization

    if m.bound_sol_pool[:cnt] >= 2 # can only warm-start the problem when pool is large enough
        ws_idx = -1
        Alp.is_min_sense(m) ? ws_obj = Inf : ws_obj = -Inf
        comp_opr = Dict(MOI.MIN_SENSE => :<, MOI.MAX_SENSE => :>)

        # Search for the pool for incumbent warm starter
        for i in 1:m.bound_sol_pool[:cnt]
            m.bound_sol_pool[:stat][i] == :Warmstarter &&
                (m.bound_sol_pool[:stat][i] = :Alive)   # reset the status if not dead
            if m.bound_sol_pool[:stat][i] != :Dead &&
               eval(comp_opr[m.sense_orig])(m.bound_sol_pool[:obj][i], ws_obj)
                ws_idx = i
                ws_obj = m.bound_sol_pool[:obj][i]
            end
        end

        if ws_idx > 0 # If a warm starter is found
            for v in m.bound_sol_pool[:vars]
                partition_cnt = length(d[v]) - 1
                active_j =
                    Alp.get_active_partition_idx(d, m.bound_sol_pool[:sol][ws_idx][v], v)
                for j in 1:partition_cnt
                    j == active_j ? JuMP.set_start_value(α[v][j], 1.0) :
                    JuMP.set_start_value(α[v][j], 0.0)
                end
            end
            m.bound_sol_pool[:stat][ws_idx] = :Warmstarter
            Alp.get_option(m, :log_level) > 0 && println(
                "!! WARM START bounding MIP using POOL SOL $(ws_idx) OBJ=$(m.bound_sol_pool[:obj][ws_idx])",
            )
        end
    end

    return
end

"""
    Method for general multilinear terms
"""
function amp_post_convhull_constrs(
    m::Optimizer,
    λ::Dict,
    α::Dict,
    indices::Any,
    dim::Tuple,
    ext_cnt::Int,
    d::Dict,
)

    # Adding λ constraints
    JuMP.@constraint(m.model_mip, sum(λ[indices][:vars]) == 1)
    JuMP.@constraint(
        m.model_mip,
        _index_to_variable_ref(m.model_mip, λ[indices][:lifted_var_idx]) ==
        dot(λ[indices][:vars], reshape(λ[indices][:vals], ext_cnt))
    )

    # Add links on each dimension
    for (cnt, i) in enumerate(indices)
        l_cnt = length(d[i])
        if m.var_type[i] == :Cont
            Alp.amp_post_inequalities_cont(m, d, λ, α, indices, dim, i, cnt)        # Add links between λ and α
        else
            error("EXCEPTION: unexpected variable type during integer related realxation")
        end
        sliced_indices =
            [collect_indices(λ[indices][:indices], cnt, [k], dim) for k in 1:l_cnt] # Add x = f(λ) for convex representation of x value
        JuMP.@constraint(
            m.model_mip,
            _index_to_variable_ref(m.model_mip, i) == sum(
                dot(
                    repeat([d[i][k]], length(sliced_indices[k])),
                    λ[indices][:vars][sliced_indices[k]],
                ) for k in 1:l_cnt
            )
        )
    end

    return
end

"""
    Constraints for univariate quadratic terms
"""
function amp_post_convhull_constrs(
    m::Optimizer,
    λ::Dict,
    α::Dict,
    monomial_idx::Int,
    discretization::Dict,
)
    partition_cnt = length(discretization[monomial_idx]) - 1
    lambda_cnt = length(discretization[monomial_idx])

    # Adding λ constraints
    JuMP.@constraint(m.model_mip, sum(λ[monomial_idx][:vars]) == 1)
    JuMP.@constraint(
        m.model_mip,
        _index_to_variable_ref(m.model_mip, λ[monomial_idx][:lifted_var_idx]) <=
        dot(λ[monomial_idx][:vars], λ[monomial_idx][:vals])
    )
    JuMP.@constraint(
        m.model_mip,
        _index_to_variable_ref(m.model_mip, λ[monomial_idx][:lifted_var_idx]) >=
        _index_to_variable_ref(m.model_mip, monomial_idx)^2
    )

    # Add SOS-2 Constraints with basic encoding
    if Alp.get_option(m, :convhull_ebd) && partition_cnt > 2
        ebd_map = embedding_map(
            lambda_cnt,
            Alp.get_option(m, :convhull_ebd_encode),
            Alp.get_option(m, :convhull_ebd_ibs),
        )
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[monomial_idx])
        for i in 1:YCnt
            JuMP.@constraint(
                m.model_mip,
                sum(λ[monomial_idx][:vars][collect(ebd_map[i])]) <= α[monomial_idx][i]
            )
            JuMP.@constraint(
                m.model_mip,
                sum(λ[monomial_idx][:vars][collect(ebd_map[i+YCnt])]) <=
                1 - α[monomial_idx][i]
            )
        end
    else
        for i in 1:lambda_cnt
            if i == 1
                JuMP.@constraint(
                    m.model_mip,
                    λ[monomial_idx][:vars][i] <= α[monomial_idx][i]
                )
            elseif i == lambda_cnt
                JuMP.@constraint(
                    m.model_mip,
                    λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1]
                )
            else
                JuMP.@constraint(
                    m.model_mip,
                    λ[monomial_idx][:vars][i] <=
                    α[monomial_idx][i-1] + α[monomial_idx][i]
                )
            end
        end
        # Add x = f(α) for regulating the domains
        JuMP.@constraint(
            m.model_mip,
            _index_to_variable_ref(m.model_mip, monomial_idx) >= sum(
                α[monomial_idx][j] * discretization[monomial_idx][j] for
                j in 1:lambda_cnt-1
            )
        )
        JuMP.@constraint(
            m.model_mip,
            _index_to_variable_ref(m.model_mip, monomial_idx) <= sum(
                α[monomial_idx][j-1] * discretization[monomial_idx][j] for
                j in 2:lambda_cnt
            )
        )
    end

    # Add x = f(λ) for convex representation
    JuMP.@constraint(
        m.model_mip,
        _index_to_variable_ref(m.model_mip, monomial_idx) ==
        dot(λ[monomial_idx][:vars], discretization[monomial_idx])
    )

    return
end

"""
    Method for multilinear terms with only continuous variables
"""
function amp_post_inequalities_cont(
    m::Optimizer,
    discretization::Dict,
    λ::Dict,
    α::Dict,
    ml_indices::Any,
    dim::Tuple,
    var_ind::Int,
    cnt::Int,
)
    lambda_cnt = length(discretization[var_ind])
    partition_cnt = lambda_cnt - 1

    # Embedding formulation
    if Alp.get_option(m, :convhull_formulation) == "sos2" &&
       Alp.get_option(m, :convhull_ebd) &&
       partition_cnt > 2
        ebd_map = embedding_map(
            lambda_cnt,
            Alp.get_option(m, :convhull_ebd_encode),
            Alp.get_option(m, :convhull_ebd_ibs),
        )
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[var_ind])
        for i in 1:YCnt
            p_sliced_indices = Alp.collect_indices(
                λ[ml_indices][:indices],
                cnt,
                collect(ebd_map[i]),
                dim,
            )
            n_sliced_indices = Alp.collect_indices(
                λ[ml_indices][:indices],
                cnt,
                collect(ebd_map[i+YCnt]),
                dim,
            )
            JuMP.@constraint(
                m.model_mip,
                sum(λ[ml_indices][:vars][p_sliced_indices]) <= α[var_ind][i]
            )
            JuMP.@constraint(
                m.model_mip,
                sum(λ[ml_indices][:vars][n_sliced_indices]) <= 1 - α[var_ind][i]
            )
        end
        Alp.get_option(m, :convhull_ebd_link) && Alp.ebd_link_xα(
            m,
            α[var_ind],
            lambda_cnt,
            discretization[var_ind],
            ebd_map[:H_orig],
            var_ind,
        )
        return
    end

    # SOS-2 Formulation
    if Alp.get_option(m, :convhull_formulation) == "sos2"
        for j in 1:lambda_cnt
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j], dim)
            if (j == 1)
                JuMP.@constraint(
                    m.model_mip,
                    sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][j]
                )
            elseif (j == lambda_cnt)
                JuMP.@constraint(
                    m.model_mip,
                    sum(λ[ml_indices][:vars][sliced_indices]) <=
                    α[var_ind][partition_cnt]
                )
            else
                JuMP.@constraint(
                    m.model_mip,
                    sum(λ[ml_indices][:vars][sliced_indices]) <= sum(α[var_ind][(j-1):j])
                )
            end
        end
        return
    elseif Alp.get_option(m, :convhull_formulation) == "facet"
        for j in 1:(partition_cnt-1) # Constraint cluster of α >= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
            JuMP.@constraint(
                m.model_mip,
                sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices])
            )
        end
        for j in 1:(partition_cnt-1) # Constraint cluster of α <= f(λ)
            sliced_indices =
                collect_indices(λ[ml_indices][:indices], cnt, [1:(j+1);], dim)
            JuMP.@constraint(
                m.model_mip,
                sum(α[var_ind][1:j]) <= sum(λ[ml_indices][:vars][sliced_indices])
            )
        end
        return
    else
        error("Must indicate a choice of convex hull formulation: sos2, facet")
    end

    return
end

function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Array, dim::Tuple)
    k = 0
    indices =
        Vector{Int}(undef, Int(prod(dim) / dim[fixed_dim] * length(fixed_partition)))
    for i in 1:prod(dim)
        ind = Tuple(CartesianIndices(l)[i])
        if ind[fixed_dim] in fixed_partition
            k += 1
            indices[k] = i
        end
    end

    return indices
end

"""
    _add_multilinear_linking_constraints(m::Optimizer, λ::Dict)

This internal function adds linking constraints between λ multipliers corresponding to multilinear terms
that share more than two variables and are partitioned. For example, suppose we have λ[i], λ[j], and 
λ[k] where i=(1,2,3), j=(1,2,4), and k=(1,2,5). λ[i] contains all multipliers for the extreme points 
in the space of (x1,x2,x3). λ[j] contains all multipliers for the extreme points in the space of (x1,x2,x4).
λ[k] contains all multipliers for the extreme points in the space of (x1,x2,x5).

Using λ[i], λ[j], or λ[k], we can express multilinear function x1*x2.
We define a linking variable μ(1,2) that represents the value of x1*x2.
Linking constraints are
    μ(1,2) == convex combination expr for x1*x2 using λ[i],    
    μ(1,2) == convex combination expr for x1*x2 using λ[j], and   
    μ(1,2) == convex combination expr for x1*x2 using λ[k].    

Thus, these constraints link between λ[i], λ[j], and λ[k] variables.

Reference: J. Kim, J.P. Richard, M. Tawarmalani, Piecewise Polyhedral Relaxations of Multilinear Optimization, 
http://www.optimization-online.org/DB_HTML/2022/07/8974.html
"""
function _add_multilinear_linking_constraints(m::Optimizer, λ::Dict)
    if isnothing(m.linking_constraints_info)
        m.linking_constraints_info = Alp._get_shared_multilinear_terms_info(
            λ,
            m.options.linking_constraints_degree_limit,
        )
    end

    if isnothing(m.linking_constraints_info)
        return
    end

    # Additional linking variables (μ) to the MIP model to keep the constraints sparse
    linking_variables =
        JuMP.@variable(m.model_mip, [i in keys(m.linking_constraints_info)])

    # Add linking constraints to the MIP model
    for (shared_multilinear_idx, multilinear_terms_idx) in m.linking_constraints_info,
        (i, multilinear_idx) in enumerate(multilinear_terms_idx)

        var_location = [findfirst(multilinear_idx .== i) for i in shared_multilinear_idx]

        lambda_expr = 0
        for idx in Iterators.product([1:d for d in λ[multilinear_idx][:dim]]...)
            var = λ[multilinear_idx][:vars][λ[multilinear_idx][:indices][idx...]]
            partitions_info = prod(
                m.discretization[i][idx[j]] for
                (i, j) in zip(shared_multilinear_idx, var_location)
            )
            lambda_expr += partitions_info * var
        end

        JuMP.@constraint(
            m.model_mip,
            linking_variables[shared_multilinear_idx] == lambda_expr
        )
    end
end

"""
    _get_shared_multilinear_terms_info(λ, linking_constraints_degree_limit)

This function checks to see if linking constraints are 
necessary for a given vector of each multilinear terms and returns the approapriate 
linking constraints information.
"""
function _get_shared_multilinear_terms_info(
    λ::Dict,
    linking_constraints_degree_limit::Union{Nothing,T} where {T<:Int64} = nothing,
)

    # Compute maximum degree of multilinear terms and return if bilinear
    max_degree = maximum([length(k) for k in keys(λ)])

    if max_degree <= 2
        return (linking_constraints_info = nothing)
    end

    # Limit the linking constraints to a prescribed multilinear degree 
    if !isnothing(linking_constraints_degree_limit) &&
       (linking_constraints_degree_limit < max_degree)
        max_degree = linking_constraints_degree_limit
    end

    # Collect all variable indices appearing in multilinear terms
    all_variables_idx = collect(union((keys(λ) |> collect)...)) |> sort

    linking_constraints_info = Dict(
        shared_multilinear_idx =>
            filter(r -> issubset(shared_multilinear_idx, r), keys(λ)) for
        deg in 2:(max_degree-1) for
        shared_multilinear_idx in Combinatorics.combinations(all_variables_idx, deg)
    )

    filter!(r -> length(r.second) >= 2, linking_constraints_info)
    if isempty(linking_constraints_info)
        return (linking_constraints_info = nothing)
    end

    return (linking_constraints_info = linking_constraints_info)
end
