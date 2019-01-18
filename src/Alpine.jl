__precompile__()

module Alpine

PODDEBUG = false

using JuMP
using MathProgBase



if VERSION < v"0.7.0-"
    import Compat: undef
    import Compat: @warn
    import Compat: occursin
    import Compat: printstyled
    import Compat: round

end

if VERSION > v"0.7.0-"
    using LinearAlgebra
    using Random
    using Distributed
    using Statistics
    using Printf
    using SparseArrays

    function ind2sub(x, i)
        return Tuple(CartesianIndices(x)[i])
    end

end




include("const.jl")

# Engine for High-level Algorithmic Control and User-interface
include("solver.jl")
include("mpb2moi.jl") # Transition file

# Engine for expression handling
include("nlexpr.jl")
include("operators.jl")

# Main Algorithmic Process
include("algorithm.jl")
include("presolve.jl")
include("amp.jl")
include("embedding.jl")
include("heuristics.jl")

# Convexification method
include("multi.jl")
include("tmc.jl")

# Model Manipulation and utilities
include("bounds.jl")
include("utility.jl")

# Othes
include("log.jl")

end # module
