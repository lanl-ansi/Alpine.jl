#=====================----- Developer's note's -----=====================#
#  Code used to use convexhull representation to convexify the non-convex model.
#========================================================================#
function amp_post_convhull(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables

    # Construct λ variable space
    for bi in keys(m.nonlinear_info)
        nl_type = m.nonlinear_info[bi][:nonlinear_type]
        if ((nl_type == :multilinear) || (nl_type == :bilinear) || (nl_type == :monomial)) && (m.nonlinear_info[bi][:convexified] == false)
            m.nonlinear_info[bi][:convexified] = true  # Bookeeping the examined terms
            ml_indices, dim, extreme_point_cnt = amp_convhull_prepare(discretization, bi)   # convert key to easy read mode
            λ = amp_convhull_λ(m, bi, ml_indices, λ, extreme_point_cnt, dim)
            λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, ones(Int,length(dim)))
            α = amp_convhull_α(m, ml_indices, α, dim, discretization)
            amp_post_convhull_constrs(m, λ, α, ml_indices, dim, extreme_point_cnt, discretization)
        end
    end

    return
end

"""
    TODO: docstring
    This function is very important.
"""
function amp_convhull_prepare(discretization::Dict, nonlinear_key::Any)

    id = Set()                      # Coverting the nonlinear indices into a different space
    for var in nonlinear_key        # This output regulates the sequence of how composing variable should be arranged
        push!(id, var.args[2])
    end

    if length(id) < length(nonlinear_key) # Got repeating terms, now the sequence matters
        id = []
        for var in nonlinear_key
            push!(id, var.args[2])
        end
    end

    dim = []
    for i in id                     # Critical!!! Ensure the same sequence
        push!(dim, length(discretization[i]))
    end


    return id, tuple([i for i in dim]...), prod(dim)
end

function amp_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, ml_indices::Any, λ::Dict, extreme_point_cnt::Int, dim::Tuple; kwargs...)

    lifted_var_idx = m.nonlinear_info[nonlinear_key][:lifted_var_ref].args[2]
    @assert !(lifted_var_idx in keys(λ))
    λ[ml_indices] = Dict(:dim=>dim,
                         :lifted_var_idx=>lifted_var_idx,
                         :indices=>reshape([1:extreme_point_cnt;], dim),
                         :vars=>@variable(m.model_mip, [1:extreme_point_cnt], lowerbound=0, upperbound=1, basename="L$(lifted_var_idx)"),
                         :vals=>ones(dim))

    return λ
end

"""
    TODO: docstring
"""
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

"""
    Less memeory & time efficient, but a easier implementation
"""
function _populate_convhull_extreme_values(discretization::Dict, ml_indices::Any, λ::Dict, extreme_point_cnt::Int)

    var_indices = collect(ml_indices)
    for i in 1:extreme_point_cnt
        sub = ind2sub(λ[ml_indices][:indices], i)
        k = 0
        for var in ml_indices
            k += 1
            @assert var == var_indices[k]
            λ[ml_indices][:vals][CartesianIndex(sub)] *= discretization[var][sub[k]]
        end
    end

    return λ
end

function amp_convhull_α(m::PODNonlinearModel, ml_indices::Any, α::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in ml_indices
        if !(i in keys(α))
            partition_cnt = length(discretization[i]) - 1
            α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(i)")
            @constraint(m.model_mip, sum(α[i]) == 1)
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

        partition_cnt = length(α[i])
        lambda_cnt = length(discretization[i])
        @assert lambda_cnt == partition_cnt + 1

        # Add links between λ and α
        valid_inequalities(m, discretization, λ, α, ml_indices, dim, i, cnt)

        # Add x = f(λ) for convex representation of x value
        sliced_indices = [collect_indices(λ[ml_indices][:indices], cnt, [k], dim) for k in 1:lambda_cnt]
        @constraint(m.model_mip, Variable(m.model_mip, i) == sum(dot(repmat([discretization[i][k]],length(sliced_indices[k])), λ[ml_indices][:vars][sliced_indices[k]]) for k in 1:lambda_cnt))
    end

    return
end

# Valid inequalities proposed when Jeff L. was here
function valid_inequalities(m::PODNonlinearModel, discretization::Dict, λ::Dict, α::Dict, ml_indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    partition_cnt = length(α[var_ind])
    lambda_cnt = length(discretization[var_ind])

    # Constriant cluster of α <= f(λ)
    for j in 1:partition_cnt
        for i in 1:max(1, min(partition_cnt-j+1, m.convhull_sweep_limit)) # At least one
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j:(j+i);], dim)
            @constraint(m.model_mip, sum(α[var_ind][j:(j+i-1)]) <= sum(λ[ml_indices][:vars][sliced_indices]))
        end
    end

    # Constraint cluster of α >= f(λ)
    for j in 1:min(partition_cnt, m.convhull_sweep_limit) # Construct cuts by sweeping in both directions
        sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
        @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [(lambda_cnt-j+1):(lambda_cnt);], dim)
        @constraint(m.model_mip, sum(α[var_ind][(dim[cnt]-j):(dim[cnt]-1)]) >= sum(λ[ml_indices][:vars][sliced_indices]))
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
