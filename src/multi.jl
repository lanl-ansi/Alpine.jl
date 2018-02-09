#=====================----- Developer's note's -----=====================#
#  Code used to use convexhull representation to convexify the non-convex model.
#========================================================================#
function amp_post_convhull(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables
    β = Dict()  # Lifted variables for exact formulation

    # Construct λ variable space
    for k in keys(m.nonlinear_terms)
        nl_type = m.nonlinear_terms[k][:nonlinear_type]
        if ((nl_type == :MULTILINEAR) || (nl_type == :BILINEAR)) && (m.nonlinear_terms[k][:convexified] == false)
            λ, α = amp_convexify_multilinear(m, k, λ, α, discretization)
        elseif nl_type == :MONOMIAL && !m.nonlinear_terms[k][:convexified]
            λ, α = amp_convexify_monomial(m, k, λ, α, discretization)
        elseif nl_type == :BINLIN && !m.nonlinear_terms[k][:convexified]
            β = amp_convexify_binlin(m, k, β)
        elseif nl_type == :BINPROD && !m.nonlinear_terms[k][:convexified]
            β = amp_convexify_binprod(m, k, β)
        elseif nl_type == :INTPROD && !m.nonlinear_terms[k][:convexified]
            λ, α = amp_convexify_intprod(m, k, λ, α, discretization)
        elseif nl_type == :BININT && !m.nonlinear_terms[k][:convexified]
            β = amp_convexify_binint(m, k, β)
        elseif nl_type == :INTLIN && !m.nonlinear_terms[k][:convexified]
            λ, α = amp_convexify_intlin(m, k, λ, α, discretization)
        elseif nl_type in [:sin, :cos] && !m.nonlinear_terms[l][:convexified]
            β = amp_convexify_sincos(m, k, β)
        end
    end

    return
end

function amp_convexify_intprod(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, discretization::Dict)

    error("No method implemented to relax INTPROD term")

    return λ, α
end

function amp_convexify_binint(m::PODNonlinearModel, k::Any, β::Dict)

    m.nonlinear_terms[k][:convexified] = true  # Bookeeping the convexified terms

    @assert length(m.nonlinear_terms[k][:var_idxs]) == 2

    lift_idx = m.nonlinear_terms[k][:y_idx]

    if haskey(β, lift_idx)
        return β
    else
        β[lift_idx] = Variable(m.model_mip, lift_idx)
    end

    bin_idx = [i for i in m.nonlinear_terms[k][:var_idxs] if m.var_type[i] == :Bin]
    int_idx = [i for i in m.nonlinear_terms[k][:var_idxs] if m.var_type[i] == :Int]

    @assert length(bin_idx) == length(int_idx) == 1

    bin_idx = bin_idx[1]
    cont_idx = cont_idx[1]

    mccormick_binlin(m.model_mip, Variable(m.model_mip, lift_idx),
        Variable(m.model_mip, bin_idx), Variable(m.model_mip, cont_idx),
        m.l_var_tight[cont_idx], m.u_var_tight[cont_idx])

    return β
end

function amp_convexify_intlin(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, discretization::Dict)

    error("No method implemented to relax INTLIN term")

    return λ, α
end

function amp_convexify_sincos(m::PODNonlinearModel, k::Any, β::Dict)

    lift_idx = m.nonlinear_terms[k][:y_idx]
    if haskey(β, lift_idx)
        return β
    else
        error("No method implemented to relax INTPROD term")
    end

    return β
end

function amp_convexify_multilinear(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, discretization::Dict)

    m.nonlinear_terms[k][:convexified] = true  # Bookeeping the convexified terms

    ml_indices, dim, extreme_point_cnt = amp_convhull_prepare(m, discretization, k)   # convert key to easy read mode
    λ = amp_convhull_λ(m, k, ml_indices, λ, extreme_point_cnt, dim)
    λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, ones(Int,length(dim)))
    α = amp_convhull_α(m, ml_indices, α, dim, discretization)
    amp_post_convhull_constrs(m, λ, α, ml_indices, dim, extreme_point_cnt, discretization)

    return λ, α
end

function amp_convexify_monomial(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, discretization::Dict)

    m.nonlinear_terms[k][:convexified] = true  # Bookeeping the convexified terms

    monomial_index, dim, extreme_point_cnt = amp_convhull_prepare(m, discretization, k, monomial=true)
    λ = amp_convhull_λ(m, k, monomial_index, λ, extreme_point_cnt, dim)
    λ = populate_convhull_extreme_values(m, discretization, monomial_index, λ)
    α = amp_convhull_α(m, [monomial_index], α, dim, discretization)
    amp_post_convhull_constrs(m, λ, α, monomial_index, dim, discretization)

    return λ, α
end

function amp_convexify_binlin(m::PODNonlinearModel, k::Any, β::Dict)

    m.nonlinear_terms[k][:convexified] = true  # Bookeeping the convexified terms

    @assert length(m.nonlinear_terms[k][:var_idxs]) == 2

    lift_idx = m.nonlinear_terms[k][:y_idx]

    if haskey(β, lift_idx)
        return β
    else
        β[lift_idx] = Variable(m.model_mip, lift_idx)
    end

    bin_idx = [i for i in m.nonlinear_terms[k][:var_idxs] if m.var_type[i] == :Bin]
    cont_idx = [i for i in m.nonlinear_terms[k][:var_idxs] if m.var_type[i] == :Cont]

    @assert length(bin_idx) == length(cont_idx) == 1

    bin_idx = bin_idx[1]
    cont_idx = cont_idx[1]

    mccormick_binlin(m.model_mip, Variable(m.model_mip, lift_idx),
        Variable(m.model_mip, bin_idx), Variable(m.model_mip, cont_idx),
        m.l_var_tight[cont_idx], m.u_var_tight[cont_idx])

    return β
end

function amp_convexify_binprod(m::PODNonlinearModel, k::Any, β::Dict)

    m.nonlinear_terms[k][:convexified] = true  # Bookeeping the convexified terms

    lift_idx = m.nonlinear_terms[k][:y_idx]
    if haskey(β, lift_idx)
        return β    # Already constructed
    else
        β[lift_idx] = Variable(m.model_mip, lift_idx)
    end

    z = Variable(m.model_mip, m.nonlinear_terms[k][:y_idx])
    x = [Variable(m.model_mip, i) for i in m.nonlinear_terms[k][:var_idxs]]
    for i in x
        @constraint(m.model_mip, z <= i)
    end
    @constraint(m.model_mip, z >= sum(x) - (length(x)-1))

    return β
end

function amp_convhull_prepare(m::PODNonlinearModel, discretization::Dict, nonlinear_key::Any; monomial=false)

    counted_var = []
    id = Set()                      # Coverting the nonlinear indices into a set
    for var in nonlinear_key        # This output regulates the sequence of how composing variable should be arranged
        m.var_type[var.args[2]] == :Cont && push!(id, var.args[2])
        m.var_type[var.args[2]] == :Cont && push!(counted_var, var.args[2])
    end

    if length(id) < length(counted_var) # Got repeating terms, now the sequence matters
        id = []
        for var in nonlinear_key
            m.var_type[var.args[2]] == :Cont && push!(id, var.args[2])
        end
    end

    dim = []
    for i in id
        push!(dim, length(discretization[i]))
    end

    monomial && return id[1], tuple(dim[1]), dim[1]   # One less dimension is required

    return id, tuple([i for i in dim]...), prod(dim)
end

function amp_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, ml_indices::Any, λ::Dict, extreme_point_cnt::Int, dim::Tuple; kwargs...)

    lifted_var_idx = m.nonlinear_terms[nonlinear_key][:lifted_var_ref].args[2]
    @assert !(lifted_var_idx in keys(λ))
    λ[ml_indices] = Dict(:dim=>dim,
                         :lifted_var_idx=>lifted_var_idx,
                         :indices=>reshape([1:extreme_point_cnt;], dim),
                         :vars=>@variable(m.model_mip, [1:extreme_point_cnt], lowerbound=0, upperbound=1, basename="L$(lifted_var_idx)"),
                         :vals=>ones(dim))

    return λ
end

function populate_convhull_extreme_values(m::PODNonlinearModel, discretization::Dict, monomial_index::Int, λ::Dict)
    λ[monomial_index][:vals] = [discretization[monomial_index][i]^2 for i in 1:length(discretization[monomial_index])]
    return λ
end

function populate_convhull_extreme_values(m::PODNonlinearModel, discretization::Dict, ml_indices::Any, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

    if level > length(dim)
        @assert length(ml_indices) == length(dim)
        @assert length(ml_indices) == length(locator)
        val = 1.0
        k = 0
        for i in ml_indices
            k += 1
            val *= discretization[i][locator[k]]                      # Calculate extreme point z-value
        end
        λ[ml_indices][:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i
            λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, locator, level+1)
        end
    end

    return λ
end

function amp_convhull_α(m::PODNonlinearModel, ml_indices::Any, α::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in ml_indices
        if !(i in keys(α))
            lambda_cnt = length(discretization[i])
            partition_cnt = length(discretization[i]) - 1
            if m.convhull_ebd && partition_cnt > 2
                αCnt = Int(ceil(log(2,partition_cnt)))
                α[i] = @variable(m.model_mip, [1:αCnt], Bin, basename=string("YL",i))
            else
                α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(i)")
                @constraint(m.model_mip, sum(α[i]) == 1)
                @constraint(m.model_mip, Variable(m.model_mip, i) >= sum(α[i][j]*discretization[i][j] for j in 1:lambda_cnt-1)) # Add x = f(α) for regulating the domains
                @constraint(m.model_mip, Variable(m.model_mip, i) <= sum(α[i][j-1]*discretization[i][j] for j in 2:lambda_cnt))
            end
        end
    end

    return α
end

function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, ml_indices::Any, dim::Tuple, extreme_point_cnt::Int, discretization::Dict)

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[ml_indices][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[ml_indices][:lifted_var_idx]) == dot(λ[ml_indices][:vars], reshape(λ[ml_indices][:vals], extreme_point_cnt)))

    cnt = 0
    for i in ml_indices
        cnt += 1
        amp_post_inequalities(m, discretization, λ, α, ml_indices, dim, i, cnt)        # Add links between λ and α
        lambda_cnt = length(discretization[i])
        sliced_indices = [collect_indices(λ[ml_indices][:indices], cnt, [k], dim) for k in 1:lambda_cnt] # Add x = f(λ) for convex representation of x value
        @constraint(m.model_mip, Variable(m.model_mip, i) == sum(dot(repmat([discretization[i][k]],length(sliced_indices[k])), λ[ml_indices][:vars][sliced_indices[k]]) for k in 1:lambda_cnt))
    end

    return
end

function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, monomial_idx::Int, dim::Tuple, discretization::Dict)

    partition_cnt = length(discretization[monomial_idx])-1
    lambda_cnt = length(discretization[monomial_idx])

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[monomial_idx][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[monomial_idx][:lifted_var_idx]) <= dot(λ[monomial_idx][:vars], λ[monomial_idx][:vals]))
    @constraint(m.model_mip, Variable(m.model_mip, λ[monomial_idx][:lifted_var_idx]) >= Variable(m.model_mip, monomial_idx)^2)

    # Add SOS-2 Constraints with basic encoding
    if m.convhull_ebd && partition_cnt > 2
        ebd_map = embedding_map(lambda_cnt, m.convhull_ebd_encode, m.convhull_ebd_ibs)
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
        @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) >= sum(α[monomial_idx][j]*discretization[monomial_idx][j] for j in 1:lambda_cnt-1))
        @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) <= sum(α[monomial_idx][j-1]*discretization[monomial_idx][j] for j in 2:lambda_cnt))
    end

    # Equalivent SOS-2 : A different encoding (solwer performance)
    # for i in 1:partition_cnt
    #     @constraint(m.model_mip, α[monomial_idx][i] <= λ[monomial_idx][:vars][i] + λ[monomial_idx][:vars][i+1])
    # end
    # @constraint(m.model_mip, α[monomial_idx][1] >= λ[monomial_idx][:vars][1])
    # @constraint(m.model_mip, α[monomial_idx][end] >= λ[monomial_idx][:vars][end])

    # Add x = f(λ) for convex representation
    @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) == dot(λ[monomial_idx][:vars], discretization[monomial_idx]))

    return
end

function amp_post_inequalities(m::PODNonlinearModel, discretization::Dict, λ::Dict, α::Dict, ml_indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    lambda_cnt = length(discretization[var_ind])
    partition_cnt = lambda_cnt - 1

    # Embedding formulation
    if m.convhull_formulation == "sos2" && m.convhull_ebd && partition_cnt > 2
        ebd_map = embedding_map(lambda_cnt, m.convhull_ebd_encode, m.convhull_ebd_ibs)
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[var_ind])
        for i in 1:YCnt
            p_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, collect(ebd_map[i]), dim)
            n_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, collect(ebd_map[i+YCnt]), dim)
            @constraint(m.model_mip, sum(λ[ml_indices][:vars][p_sliced_indices]) <= α[var_ind][i])
            @constraint(m.model_mip, sum(λ[ml_indices][:vars][n_sliced_indices]) <= 1-α[var_ind][i])
        end
        m.convhull_ebd_link && ebd_link_xα(m, α[var_ind], lambda_cnt, discretization[var_ind], ebd_map[:H_orig], var_ind)
        return
    end

    # SOS-2 Formulation
    if m.convhull_formulation == "sos2"
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
    elseif m.convhull_formulation == "facet"
        for j in 1:(partition_cnt-1) # Constraint cluster of α >= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        for j in 1:(partition_cnt-1) # Constriant cluster of α <= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:(j+1);], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) <= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        return
    elseif m.convhull_formulation == "mini"
        for j in 1:min(partition_cnt, 1) # Constraint cluster of α >= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1;], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [lambda_cnt;], dim)
            @constraint(m.model_mip, sum(α[var_ind][(dim[cnt]-j):(dim[cnt]-1)]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        for j in 1:partition_cnt         # Constriant cluster of α <= f(λ)
            for i in 1:max(1, min(partition_cnt-j+1, 1)) # At least one
                sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j:(j+i);], dim)
                @constraint(m.model_mip, sum(α[var_ind][j:(j+i-1)]) <= sum(λ[ml_indices][:vars][sliced_indices]))
            end
        end
        return
    else
        error("Must indicate a choice of convex hull formulation. ?(mini, sos2, facet)")
    end

    return
end

function collect_indices(l::Array, locator::Tuple, dim::Tuple)

    k = 0
    indices = Vector{Int}(2^length(dim))
    for i in 1:prod(dim)
        ind = ind2sub(l, i)
        diff = [((ind[i] - locator[i] == 0) || (ind[i] - locator[i] == 1)) for i in 1:length(dim)]
        if prod(diff)
            k +=1
            indices[k] = i
        end
    end

    return indices
end

function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Array, dim::Tuple)

	k = 0
	indices = Vector{Int}(Int(prod(dim)/dim[fixed_dim]*length(fixed_partition)))
	for i in 1:prod(dim)
		ind = ind2sub(l, i)
		if ind[fixed_dim] in fixed_partition
			k += 1
			indices[k] = i
		end
	end

	return indices
end

function resolve_lifted_var_value(m::PODNonlinearModel, sol_vec::Array)

    @assert length(sol_vec) == m.num_var_orig
    sol_vec = [sol_vec; fill(NaN, m.num_var_linear_mip+m.num_var_nonlinear_mip)]

    for i in 1:length(m.nonlinear_terms)
        for bi in keys(m.nonlinear_terms)
            if m.nonlinear_terms[bi][:id] == i
                lvar_idx = m.num_var_orig + i
                if haskey(m.nonlinear_terms[bi], :evaluator)
                    sol_vec[lvar_idx] = m.nonlinear_terms[bi][:evaluator](m.nonlinear_terms[bi], sol_vec)
                else
                    sol_vec[lvar_idx] = 0.5*m.discretization[lvar_idx][1] + 0.5*m.discretization[lvar_idx][end]
                end
            end
        end
    end

    return sol_vec
end
