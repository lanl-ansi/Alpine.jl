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
        if ((nl_type == :multilinear) || (nl_type == :bilinear)) && (m.nonlinear_info[bi][:convexified] == false)
            m.nonlinear_info[bi][:convexified] = true  # Bookeeping the examined terms
            ml_indices, dim, extreme_point_cnt = amp_convhull_prepare(discretization, bi)   # convert key to easy read mode
            # @show bi, ml_indices, dim, extreme_point_cnt
            λ = amp_convhull_λ(m, bi, ml_indices, λ, extreme_point_cnt, dim)
            λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, ones(Int,length(dim)))
            # @show λ[ml_indices][:vals]
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
    id = Set()                 # This may be something else, TODO: be safe in regulating the shape
    dim = []
    for var in nonlinear_key        # This output regulates the sequence of how composing variable should be arranged
        push!(id, var.args[2])
    end
    for i in id                     # Ensure the same sequence
        push!(dim, length(discretization[i]))
    end
    return id, tuple([i for i in dim]...), prod(dim)
end

function amp_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, ml_indices::Set, λ::Dict, extreme_point_cnt::Int, dim::Tuple; kwargs...)

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
function populate_convhull_extreme_values(m::PODNonlinearModel, discretization::Dict, ml_indices::Set, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

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

function amp_convhull_α(m::PODNonlinearModel, ml_indices::Set, α::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in ml_indices
        if i != keys(α)
            partition_cnt = length(discretization[i]) - 1
            α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(i)")
            @constraint(m.model_mip, sum(α[i]) == 1)
            @constraint(m.model_mip, Variable(m.model_mip, i) >= dot(discretization[i][1:(end-1)], α[i]))
            @constraint(m.model_mip, Variable(m.model_mip, i) <= dot(discretization[i][2:end], α[i]))
        end
    end

    return α
end

function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, ml_indices::Set, dim::Tuple, extreme_point_cnt::Int, discretization::Dict)

    @constraint(m.model_mip, sum(λ[ml_indices][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[ml_indices][:lifted_var_idx]) == dot(λ[ml_indices][:vars], reshape(λ[ml_indices][:vals], extreme_point_cnt)))
    # Constaint 3
    cnt = 0
    for i in ml_indices
        cnt += 1
        @assert length(α[i]) == (length(discretization[i]) - 1)
        partition_cnt = length(α[i])
        for k in 1:partition_cnt
            slice_indices = collect_indices(λ[ml_indices][:indices], cnt, k, dim)
            @constraint(m.model_mip, α[i][k] <= sum(λ[ml_indices][:vars][slice_indices]))
        end
    end

    return
end


function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Int, dim::Tuple)

	k = 0
	indices = Vector{Int}(Int(prod(dim)/dim[fixed_dim]*2))
	for i in 1:prod(dim)
		ind = ind2sub(l, i)
		(fixed_partition >= dim[fixed_dim]) && error("indices collection input error, exceeding partition counts")
		if ind[fixed_dim] in [fixed_partition, fixed_partition+1]
			k += 1
			indices[k] = i
		end
	end

	return indices
end
