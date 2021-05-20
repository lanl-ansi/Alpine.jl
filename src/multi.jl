function amp_post_convhull(m::Optimizer; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_disc) ? d = options[:use_disc] : d = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables
    β = Dict()  # Lifted variables for exact formulation

    # Convexification Treatment for Complex Non-Convex Terms
    for k in keys(m.nonconvex_terms)
        nl_type = m.nonconvex_terms[k][:nonlinear_type]
        if ((nl_type == :MULTILINEAR) || (nl_type == :BILINEAR)) && (m.nonconvex_terms[k][:convexified] == false)
            λ, α = amp_convexify_multilinear(m, k, λ, α, d)
        elseif nl_type == :MONOMIAL && !m.nonconvex_terms[k][:convexified]
            λ, α = amp_convexify_monomial(m, k, λ, α, d)
        elseif nl_type == :BINLIN && !m.nonconvex_terms[k][:convexified]
            β = amp_convexify_binlin(m, k, β)
        elseif nl_type == :BINPROD && !m.nonconvex_terms[k][:convexified]
            β = amp_convexify_binprod(m, k, β)
        elseif nl_type == :BININT && !m.nonconvex_terms[k][:convexified]
            β = amp_convexify_binint(m, k, β)
        elseif nl_type == :INTPROD && !m.nonconvex_terms[k][:convexified]
            error("Integer features are OFF. No support for INTPROD at this moment.")
        elseif nl_type == :INTLIN && !m.nonconvex_terms[k][:convexified]
            error("Integer features are OFF. No support for INTLIN at this moment.")
        elseif nl_type in [:sin, :cos] && !m.nonconvex_terms[k][:convexified]
            error("Integer features are OFF. No support for INTLIN at this moment.")
        end
    end

    # Experimental code for Warm starting
    get_option(m, :convhull_warmstart) && !isempty(m.best_bound_sol) && amp_warmstart_α(m, α)

    return
end

function amp_convexify_multilinear(m::Optimizer, k::Any, λ::Dict, α::Dict, discretization::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    ml_indices, dim, extreme_point_cnt = amp_convhull_prepare(m, discretization, k)   # convert key to easy read mode
    λ = amp_convhull_λ(m, k, ml_indices, λ, extreme_point_cnt, dim)
    λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, ones(Int,length(dim)))
    α = amp_convhull_α(m, ml_indices, α, dim, discretization)
    amp_post_convhull_constrs(m, λ, α, ml_indices, dim, extreme_point_cnt, discretization)

    return λ, α
end

function amp_convexify_monomial(m::Optimizer, k::Any, λ::Dict, α::Dict, discretization::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    monomial_index, dim, extreme_point_cnt = amp_convhull_prepare(m, discretization, k, monomial=true)
    λ = amp_convhull_λ(m, k, monomial_index, λ, extreme_point_cnt, dim)
    λ = populate_convhull_extreme_values(m, discretization, monomial_index, λ, 2)
    α = amp_convhull_α(m, [monomial_index], α, dim, discretization)
    amp_post_convhull_constrs(m, λ, α, monomial_index, dim, discretization)

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

    mccormick_binlin(m.model_mip, _index_to_variable_ref(m.model_mip, lift_idx),
        _index_to_variable_ref(m.model_mip, bin_idx), _index_to_variable_ref(m.model_mip, cont_idx),
        m.l_var_tight[cont_idx], m.u_var_tight[cont_idx])

    return β
end

amp_convexify_binint(m::Optimizer, k::Any, β::Dict) = amp_convexify_binlin(m, k, β)

function amp_convexify_binprod(m::Optimizer, k::Any, β::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    lift_idx = m.nonconvex_terms[k][:y_idx]
    if haskey(β, lift_idx)
        return β    # Already constructed
    else
        β[lift_idx] = _index_to_variable_ref(m.model_mip, lift_idx)
    end

    z = _index_to_variable_ref(m.model_mip, m.nonconvex_terms[k][:y_idx])
    x = [_index_to_variable_ref(m.model_mip, i) for i in m.nonconvex_terms[k][:var_idxs]]
    for i in x
        @constraint(m.model_mip, z <= i)
    end
    @constraint(m.model_mip, z >= sum(x) - (length(x)-1))

    return β
end

"""
    Method for general nonlinear terms
"""
function amp_convhull_prepare(m::Optimizer, d::Dict, nonlinear_key::Any; monomial=false)

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
function amp_convhull_λ(m::Optimizer, nonlinear_key::Any, indices::Any, λ::Dict, ext_cnt::Int, dim::Tuple)

    y_idx = m.nonconvex_terms[nonlinear_key][:y_idx]

    @assert !(y_idx in keys(λ))
    λ[indices] = Dict(:dim=>dim,
                     :lifted_var_idx=>y_idx,
                     :indices=>reshape([1:ext_cnt;], dim),
                     :vars=>@variable(m.model_mip, [1:ext_cnt], lower_bound=0, base_name="L$(y_idx)"),
                     :vals=>ones(dim))

    return λ
end

"""
    Method for power terms
"""
function populate_convhull_extreme_values(m::Optimizer, d::Dict, mono_idx::Int, λ::Dict, p::Int)
    λ[mono_idx][:vals] = [d[mono_idx][i]^p for i in 1:length(d[mono_idx])]
    return λ
end

"""
    Method for regular muiltilinear terms
"""
function populate_convhull_extreme_values(m::Optimizer, discretization::Dict, indices::Any, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

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
            λ = populate_convhull_extreme_values(m, discretization, indices, λ, dim, locator, level+1)
        end
    end

    return λ
end

"""
    General Method for all term
"""
function amp_convhull_α(m::Optimizer, indices::Any, α::Dict, dim::Tuple, discretization::Dict)

    for i in indices
        if !(i in keys(α))
            lambda_cnt = length(discretization[i])
            partition_cnt = length(discretization[i]) - 1
            if get_option(m, :convhull_ebd) && partition_cnt > 2
                αCnt = Int(ceil(log(2,partition_cnt)))
                α[i] = @variable(m.model_mip, [1:αCnt], Bin, base_name=string("YL",i))
            else
                α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, base_name="A$(i)")
                @constraint(m.model_mip, sum(α[i]) == 1)
                @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, i) >= sum(α[i][j]*discretization[i][j] for j in 1:lambda_cnt-1)) # Add x = f(α) for regulating the domains
                @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, i) <= sum(α[i][j-1]*discretization[i][j] for j in 2:lambda_cnt))
            end
        end
    end

    return α
end

amp_convhull_α(m::Optimizer, idx::Int, α::Dict, dim, d::Dict) = amp_convhull_α(m, [idx], α, dim, d)

function amp_no_good_cut_α(m::Optimizer, α::Dict)

    println("Global Incumbent solution objective = $(m.best_obj)")

    for i in 1:m.bound_sol_pool[:cnt]
        (m.bound_sol_pool[:stat][i] == :Cutoff) && (m.bound_sol_pool[:stat][i] = :Alive)
        if m.best_obj < m.bound_sol_pool[:obj][i] && m.bound_sol_pool[:stat][i] == :Alive
            no_good_idxs = keys(m.bound_sol_pool[:disc][i])
            no_good_size = length(no_good_idxs) - 1
            @constraint(m.model_mip, sum(α[v][m.bound_sol_pool[:disc][i][v]] for v in no_good_idxs) <= no_good_size)
            get_option(m, :log_level) > 0 && println("!! GLOBAL cuts off POOL_SOL-$(i) POOL_OBJ=$(m.bound_sol_pool[:obj][i])!")
            m.bound_sol_pool[:stat][i] = :Cutoff
        end
    end

    return
end

function amp_warmstart_α(m::Optimizer, α::Dict)

    d = m.discretization

    if m.bound_sol_pool[:cnt] >= 2 # can only warm-start the problem when pool is large enough
        ws_idx = -1
        is_min_sense(m) ? ws_obj = Inf : ws_obj = -Inf
        comp_opr = Dict(MOI.MIN_SENSE => :<, MOI.MAX_SENSE => :>)

        # Search for the pool for incumbent warm starter
        for i in 1:m.bound_sol_pool[:cnt]
            m.bound_sol_pool[:stat][i] == :Warmstarter && (m.bound_sol_pool[:stat][i] = :Alive)   # reset the status if not dead
            if m.bound_sol_pool[:stat][i] != :Dead && eval(comp_opr[m.sense_orig])(m.bound_sol_pool[:obj][i], ws_obj)
                ws_idx = i
                ws_obj = m.bound_sol_pool[:obj][i]
            end
        end

        if ws_idx > 0 # If a warm starter is found
            for v in m.bound_sol_pool[:vars]
                partition_cnt = length(d[v])-1
                active_j = get_active_partition_idx(d, m.bound_sol_pool[:sol][ws_idx][v], v)
                for j = 1:partition_cnt
                    j == active_j ? set_start_value(α[v][j], 1.0) : set_start_value(α[v][j], 0.0)
                end
            end
            m.bound_sol_pool[:stat][ws_idx] = :Warmstarter
            get_option(m, :log_level) > 0 && println("!! WARM START bounding MIP using POOL SOL $(ws_idx) OBJ=$(m.bound_sol_pool[:obj][ws_idx])")
        end
    end

    return
end

"""
    Method for general multilinear terms with/without integer variables
"""
function amp_post_convhull_constrs(m::Optimizer, λ::Dict, α::Dict, indices::Any, dim::Tuple, ext_cnt::Int, d::Dict)

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[indices][:vars]) == 1)
    @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, λ[indices][:lifted_var_idx]) == dot(λ[indices][:vars], reshape(λ[indices][:vals], ext_cnt)))

    # Add links on each dimension
    for (cnt, i) in enumerate(indices)
        l_cnt = length(d[i])
        if m.var_type[i] == :Cont
            amp_post_inequalities_cont(m, d, λ, α, indices, dim, i, cnt)        # Add links between λ and α
        else
            error("EXCEPTION: unexpected variable type during integer related realxation")
        end
        sliced_indices = [collect_indices(λ[indices][:indices], cnt, [k], dim) for k in 1:l_cnt] # Add x = f(λ) for convex representation of x value
        @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, i) == sum(dot(repeat([d[i][k]],length(sliced_indices[k])), λ[indices][:vars][sliced_indices[k]]) for k in 1:l_cnt))
    end

    return
end

"""
    Method for power-2 term
"""
function amp_post_convhull_constrs(m::Optimizer, λ::Dict, α::Dict, monomial_idx::Int, dim::Tuple, discretization::Dict)

    partition_cnt = length(discretization[monomial_idx])-1
    lambda_cnt = length(discretization[monomial_idx])

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[monomial_idx][:vars]) == 1)
    @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, λ[monomial_idx][:lifted_var_idx]) <= dot(λ[monomial_idx][:vars], λ[monomial_idx][:vals]))
    @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, λ[monomial_idx][:lifted_var_idx]) >= _index_to_variable_ref(m.model_mip, monomial_idx)^2)

    # Add SOS-2 Constraints with basic encoding
    if get_option(m, :convhull_ebd) && partition_cnt > 2
        ebd_map = embedding_map(lambda_cnt, get_option(m, :convhull_ebd_encode), get_option(m, :convhull_ebd_ibs))
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[monomial_idx])
        for i in 1:YCnt
            @constraint(m.model_mip, sum(λ[monomial_idx][:vars][collect(ebd_map[i])]) <= α[monomial_idx][i])
            @constraint(m.model_mip, sum(λ[monomial_idx][:vars][collect(ebd_map[i+YCnt])]) <= 1-α[monomial_idx][i])
        end
    else
        for i in 1:lambda_cnt
            if i == 1
                @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i])
            elseif i == lambda_cnt
                @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1])
            else
                @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1] + α[monomial_idx][i])
            end
        end
        # Add x = f(α) for regulating the domains
        @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, monomial_idx) >= sum(α[monomial_idx][j]*discretization[monomial_idx][j] for j in 1:lambda_cnt-1))
        @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, monomial_idx) <= sum(α[monomial_idx][j-1]*discretization[monomial_idx][j] for j in 2:lambda_cnt))
    end

    # Add x = f(λ) for convex representation
    @constraint(m.model_mip, _index_to_variable_ref(m.model_mip, monomial_idx) == dot(λ[monomial_idx][:vars], discretization[monomial_idx]))

    return
end

"""
    Method for regular multilinear terms (terms that only has continuous variables)
"""
function amp_post_inequalities_cont(m::Optimizer, discretization::Dict, λ::Dict, α::Dict, ml_indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    lambda_cnt = length(discretization[var_ind])
    partition_cnt = lambda_cnt - 1

    # Embedding formulation
    if get_option(m, :convhull_formulation) == "sos2" && get_option(m, :convhull_ebd) && partition_cnt > 2
        ebd_map = embedding_map(lambda_cnt, get_option(m, :convhull_ebd_encode), get_option(m, :convhull_ebd_ibs))
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[var_ind])
        for i in 1:YCnt
            p_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, collect(ebd_map[i]), dim)
            n_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, collect(ebd_map[i+YCnt]), dim)
            @constraint(m.model_mip, sum(λ[ml_indices][:vars][p_sliced_indices]) <= α[var_ind][i])
            @constraint(m.model_mip, sum(λ[ml_indices][:vars][n_sliced_indices]) <= 1-α[var_ind][i])
        end
        get_option(m, :convhull_ebd_link) && ebd_link_xα(m, α[var_ind], lambda_cnt, discretization[var_ind], ebd_map[:H_orig], var_ind)
        return
    end

    # SOS-2 Formulation
    if get_option(m, :convhull_formulation) == "sos2"
        for j in 1:lambda_cnt
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j], dim)
            if (j == 1)
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][j])
            elseif (j == lambda_cnt)
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][partition_cnt])
            else
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= sum(α[var_ind][(j-1):j]))
            end
        end
        return
    elseif get_option(m, :convhull_formulation) == "facet"
        for j in 1:(partition_cnt-1) # Constraint cluster of α >= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        for j in 1:(partition_cnt-1) # Constraint cluster of α <= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:(j+1);], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) <= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        return
    else
        error("Must indicate a choice of convex hull formulation. ?(sos2, facet)")
    end

    return
end

"""
    [Experimental Function]
    Method for multilinear terms with discrete variables
"""
function amp_post_inequalities_int(m::Optimizer, d::Dict, λ::Dict, α::Dict, indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    l_cnt = length(d[var_ind])
    p_cnt = l_cnt - 1

    # Embedding formulation
    get_option(m, :convhull_ebd) && @warn "Embedding is currently not supported for multilinear terms with discrete variables"

    # SOS-2 Formulation
    if get_option(m, :convhull_formulation) == "sos2"
        for j in 1:l_cnt
            sliced_indices = collect_indices(λ[indices][:indices], cnt, [j], dim)
            if (j == 1)
                rate_r = amp_pick_ratevec(d[var_ind], j)
                @constraint(m.model_mip, sum(λ[indices][:vars][sliced_indices]) <= α[var_ind][j] * rate_r)
            elseif (j == l_cnt)
                rate_l = amp_pick_ratevec(d[var_ind], j)
                @constraint(m.model_mip, sum(λ[indices][:vars][sliced_indices]) <= α[var_ind][p_cnt] * rate_l)
            else
                rate_l, rate_r = amp_pick_ratevec(d[var_ind], j)
                @constraint(m.model_mip, sum(λ[indices][:vars][sliced_indices]) <= α[var_ind][j-1] * rate_l + α[var_ind][j] * rate_r)
            end
        end
    else
        error("Only SOS-2 formulation supports convexification of terms with integer variables.")
    end

    return
end

"""
    [Experimental Function]
    This is a utility function that can be used to measure the coefficients for formulations that convexify terms with integer variables.
"""
function amp_pick_ratevec(partvec::Vector, i::Int)

    λCnt = length(partvec)

    if i == 1
        mod(partvec[i], 0.5) == 0.0 && partvec[i+1] - partvec[i] == 1.0 && return 0.5
    elseif i == λCnt
        mod(partvec[i], 0.5) == 0.0 && partvec[i] - partvec[i-1] == 1.0 && return 0.5
    else
        if mod(partvec[i], 0.5) == 0.0 && partvec[i] - partvec[i-1] == 1.0
            l = 0.5
        else
            l = 1.0
        end
        if mod(partvec[i], 0.5) == 0.0 && partvec[i+1] - partvec[i] == 1.0
            r = 0.5
        else
            r = 1.0
        end
        return l, r
    end

    return 1.0
end

function amp_collect_tight_regions(partvec::Vector)

    PCnt = length(partvec) - 1
    PCnt == 1 && return []

    tight_regions = []
    for i in 1:PCnt
        if i == 1
            (partvec[i]+1.0==partvec[i+1]) && (partvec[i+1]+1.0==partvec[i+2]) && push!(tight_regions,i)
        elseif i == PCnt
            (partvec[i]+1.0==partvec[i+1]) && (partvec[i]-1.0==partvec[i-1]) && push!(tight_regions, i)
        else
            (partvec[i]+1.0==partvec[i+1]) && (partvec[i+1]+1.0==partvec[i+2]) && (partvec[i]-1.0==partvec[i-1]) && push!(tight_regions, i)
        end
    end

    return tight_regions
end

function amp_post_λ_upperbound(m::Optimizer, λ::Dict, indices::Any, dim::Tuple, d::Dict, tregions::Vector, reg=[], level=0)

    if level == length(indices)
        isempty(tregions[level]) && return
        sliced_indices = Set(collect_indices(λ[indices][:indices], 1, [reg[1]; reg[1]+1], dim))
        for i in 2:length(reg)
            sliced_indices = intersect(sliced_indices, Set(collect_indices(λ[indices][:indices], i, [reg[i],reg[i]+1], dim)))
        end
        for i in sliced_indices
            JuMP.set_upper_bound(λ[indices][:vars][i], (1/2)^level)
        end
        return
    end

    for i in 1:length(tregions[level+1])
        push!(reg, tregions[level+1][i])
        amp_post_λ_upperbound(m, λ, indices, dim, d, tregions, reg, level+1)
        length(reg) < level && error("Something is wrong")
        length(reg) > level && pop!(reg)
    end

    return
end

function amp_post_λ_upperbound(m::Optimizer, λ::Dict, indices::Any, ub::Float64)

    for i in λ[indices][:vars] JuMP.set_upper_bound(i, ub) end

    return
end

function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Array, dim::Tuple)

	k = 0
	indices = Vector{Int}(undef, Int(prod(dim)/dim[fixed_dim]*length(fixed_partition)))
	for i in 1:prod(dim)
		ind = Tuple(CartesianIndices(l)[i])
		if ind[fixed_dim] in fixed_partition
			k += 1
			indices[k] = i
		end
	end

	return indices
end

# function collect_indices(l::Array, locator::Tuple, dim::Tuple)

#     k = 0
#     indices = Vector{Int}(2^length(dim))
#     for i in 1:prod(dim)
#         ind = Tuple(CartesianIndices(l)[i])
#         diff = [((ind[i] - locator[i] == 0) || (ind[i] - locator[i] == 1)) for i in 1:length(dim)]
#         if prod(diff)
#             k +=1
#             indices[k] = i
#         end
#     end

#     return indices
# end
