#=====================----- Developer's note's -----=====================#
#  Code used to use convexhull representation to convexify the non-convex model.
#========================================================================#
function amp_post_convhull(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables
    f = Dict()  # Extreme points function values

    # Construct λ variable space
    for bi in keys(m.nonlinear_info)
        nl_type = m.nonlinear_info[bi][:nonlinear_type]
        if ((nl_type == :multilinear) || (nl_type == :bilinear)) && (m.nonlinear_info[bi][:convexified] == false)
            m.nonlinear_info[bi][:convexified] = true  # Bookeeping the examined terms
            ml_indices, dim = amp_convhull_id_and_dim(discretization, bi)
            extreme_point_cnt = prod(dim)
            λ = post_convhull_λ(m, bi, λ, extreme_point_cnt, dim)
            λ = populate_extreme_values(m, ml_indices, λ, dim, ones(length(dim)))
            α = amp_post_convhull_α_vars(m, ml_indices, α, dim, m.discretization)
            # amp_post_convhull_constrs()
        end
    end

    return
end

function amp_convhull_id_and_dims(discretization::Dict, nonlinear_key::Any)
    id = Set()                 # This may be something else, TODO: be safe in regulating the shape
    dim = []
    for var in nonlinear_key
        push!(id, var.args[2])
        push!(dim, length(m.discretization[var.args[2]]))
    end
    return id, tuple([i for i in dims]...)
end

function amp_post_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, λ::Dict, extreme_point_cnt::Int, dim::Tuple; kwargs...)
    id = m.nonlinear_info[nonlinear_key][:id]
    λ[id] = Dict(:indices=>reshape(1:extreme_point_cnt, dim),
                 :vars=>@variable(m.model_mip, [1:extreme_point_cnt], lowerbound=0, upperbound=1, basename="L$(string(hash(id)))"),
                 :vals=>ones(dim))
    return λ
end

function populate_convhull_extreme_values(m::PODNonlinearModel, ml_indices::Set, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

    if level == length(dim)
        @assert length(ml_indices) == length(dim)
        @assert length(ml_indices) == length(locator)
        val = 1.0
        for i in ml_indices
            val *= m.discretization[i][locator[i]]                      # Calculate extreme point z-value
        end
        λ[:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i                                          # Fix this level's index
            λ = populate_convhull_extreme_values(m, ml_indices, λ, dim, locator, level+1)
        end
    end

    return λ
end

function amp_post_convhull_α_vars(m::PODNonlinearModel, ml_indices::Set, α::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in ml_indices
        partition_cnt = length(discretization[i]) - 1
        α[ml_indices] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(hash(ml_indices))")
    end

    return α
end

function amp_post_convhull_constrs(m::PODNonlinearModel, ml_indices::Set, α::Dict, dim::Tuple)
    # Build all necessary constriants according to user-options
    # Harhsa, SOS-1, SOS-2, Logrithm...

    # Currently, focus on SOS-1 implementation

    amp_post_convhull_one_partition(m, α)
    return
end

function amp_post_convhull_one_partition(m::PODNonlinearModel, α::Dict)
    # each variable must take only one domain
    for i in keys(α)
        @constraint(m, sum(α[i]) == 1)
    end
    return
end

"""
    Not finished.
"""
function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Int, dim::Tuple, ind::Int=1, locator=nothing, indices=[])

    if ind > length(dim)
        return indices
    else
        for i in dim(ind)
            locator = locator[ind] + 1
            if abs(i-fixed_partition) <=1 && locator[fixed_dim] == ind
                push!(indices, l[CartesianIndex(locator)])
            end
            indices = collect_indices(l, fixed_dim, fixed_partition, dim, ind+1, indices)
            locator = locator[ind] - 1
        end
    end

    return indices
end
