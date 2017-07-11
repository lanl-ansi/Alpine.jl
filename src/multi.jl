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
            λ = post_convhull_λ(m, bi, ml_indices, λ, extreme_point_cnt, dim)
            λ = populate_extreme_values(m, ml_indices, λ, dim, ones(length(dim)))
            α = amp_post_convhull_α_vars(m, ml_indices, α, dim, m.discretization)
            # amp_post_convhull_constrs()
        end
    end

    return
end

function amp_convhull_prepare(discretization::Dict, nonlinear_key::Any)
    id = Set()                 # This may be something else, TODO: be safe in regulating the shape
    dim = []
    for var in nonlinear_key
        push!(id, var.args[2])
        push!(dim, length(m.discretization[var.args[2]]))
    end
    return id, tuple([i for i in dims]...), prod(dim)
end

function amp_post_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, ml_indices::Set, λ::Dict, extreme_point_cnt::Int, dim::Tuple; kwargs...)

    lifted_var_idx = m.nonlinear_info[nonlinear_key][:id]
    @assert !(lifted_var_idx in keys(λ))  # This
    λ[lifted_var_idx] = Dict(:dim=>dim,
                             :var_indices=>ml_indices,
                             :indices=>reshape([1:extreme_point_cnt;], dim),
                             :vars=>@variable(m.model_mip, [1:extreme_point_cnt], lowerbound=0, upperbound=1, basename="L$(lifted_var_idx)"),
                             :vals=>ones(extreme_point_cnt))

    return λ
end

"""

"""
function populate_convhull_extreme_values(m::PODNonlinearModel, ml_indices::Set, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

    if level == length(dim)
        @assert length(ml_indices) == length(dim)
        @assert length(ml_indices) == length(locator)
        val = 1.0
        for i in ml_indices
            val *= m.discretization[i][locator[i]]                      # Calculate extreme point z-value
        end
        λ[id][:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i                                          # Fix this level's index
            λ = populate_convhull_extreme_values(m, ml_indices, λ, dim, locator, level+1)
        end
    end

    return reshape(λ[id][:vals], prod(dim))
end

function amp_post_convhull_α_vars(m::PODNonlinearModel, ml_indices::Set, α::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in ml_indices
        if i != keys(α)
            partition_cnt = length(discretization[i]) - 1
            α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(i))")
        end
    end

    return α
end

function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict)
    # Build all necessary constriants for convex-hull representation

    # Constraint Group A
    for i in keys(α)
        @constraint(m.model_mip, sum(α[i]) == i)
        @constraint(m, Variable(m, i) >= dot(m.discretization[i][1:(end-1)], α[i]))
        @constraint(m, Variable(m, i) <= dot(m.discretization[i][2:end], α[i]))
    end

    # Constraint Group B
    for i in keys(λ)
        @constraint(m.model_mip, sum(λ[i][:vars]) == 1)
        @constraint(m.model_mip, Variable(m.model_mip, i) == dot(λ[i][:vars], λ[i][:vals]))
        # Constaint 3
        for j in λ[i][:var_indices]
            @assert length(α[j]) == (m.discretization[j] - 1)
            partition_cnt = length(α[j])
            @constarint(m.model_mip, [k in 1:partition_cnt], α[j] <= sum(λ[collect_indices(λ[i][:indices], j, k, λ[i][:dim])]))
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
