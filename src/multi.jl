#=====================----- Developer's note's -----=====================#
#  Code used to use convexhull representation to convexify the non-convex model.
#========================================================================#
function amp_post_convhull(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points convex multiplier
    α = Dict()  # Partitioning Variables
    f = Dict()  # Extreme points function values

    # Construct λ variable space
    for bi in keys(m.nonlinear_info)
        # Downward compatability: convex hull representation is suitable for :monomial and :bilinear
        if (m.nonlinear_info[bi][:nonlinear_type] == :multilinear) || m.convex_disable_tmc * (m.nonlinear_info[bi][:nonlinear_info] in [:monomial, :bilinear])
            m.nonlinear_info[bi][:convexified] = true  # Bookeeping the examined terms
            id, dim = generate_dims(discretization, bi)
            extreme_point_cnt = prod(dim)
            λ = post_convhull_λ(m, bi[:id], λ, extreme_point_cnt, dim)
            λ = populate_extreme_values(m, id, λ, dim, ones(length(dim)))
            α = amp_post_convhull_α_vars(m, id, α, dim, m.discretization)
            # amp_post_convhull_constrs()
        end
    end

    return
end

function amp_convhull_id_and_dims(discretization::Dict, nonlinear_key::Any)
    id = Set()                  # This may be something else, TODO: be safe in regulating the shape
    dim = []
    for var in nonlinear_key
        push!(id, var.args[2])
        push!(dim, length(m.discretization[var.args[2]]))
    end
    return id, tuple([i for i in dims]...)
end

function amp_post_convhull_λ(m::PODNonlinearModel, id::Set, λ::Dict, dim::Tuple; kwargs...)
    λ[id] = Dict(:indices=>reshape(1:prod(dim), dim),
                 :vars=>@variable(m.model_mip, [1:prod(dim)], lowerbound=0, upperbound=1, basename="L$(string(hash(id)))"),
                 :vals=>ones(dim))
    return λ
end

function populate_convhull_extreme_values(m::PODNonlinearModel, id::Set, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

    if level == length(dim)
        @assert length(id) == length(dim)
        @assert length(id) == length(locator)
        val = 1.0
        for i in id
            val *= m.discretization[i][locator[i]]                      # Calculate extreme point z-value
        end
        λ[:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i                                          # Fix this level's index
            λ = populate_extreme_values(m, id, λ, dim, locator, level+1)
        end
    end

    return λ
end

function amp_post_convhull_α_vars(m::PODNonlinearModel, id::Set, α::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in id
        partition_cnt = length(discretization[i]) - 1
        α[id] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(hash(id))")
    end

    return α
end

function amp_post_convhull_constrs(m::PODNonlinearModel, id::Set, α::Dict, dim::Tuple)
    # Build all necessary constriants according to user-options
    # Harhsa, SOS-1, SOS-2, Logrithm...

    # Currently, focus on SOS-1 implementation

    # Post SOS-1
    for i in id
        # amp_post_convhull_SOS1()
        # amp_post_convhull_domain()
        # amp_post_convhull_domain()
        # amp_post_convhull_hull()
        # amp_post_convhull_locator()
    end

    return
end

function amp_post_convhull_sos1()

    return
end
