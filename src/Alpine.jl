__precompile__()

module Alpine

using JuMP
using LinearAlgebra
using Random
using Statistics

const ALPINE_DEBUG = false
const ALP = Alpine
const LA = LinearAlgebra

include("const.jl")

# Engine for High-level Algorithmic Control and User-interface
include("matrix_opt_interface.jl")
include("moi_function2expr.jl")
include("solver.jl")

# Engine for expression handling
include("nlexpr.jl")
include("operators.jl")

# Main Algorithm
include("algorithm.jl")
include("presolve.jl")
include("amp.jl")
include("embedding.jl")
include("heuristics.jl")

# Convexification 
include("multi.jl")
include("tmc.jl")

# Model manipulation and utilities
include("bounds.jl")
include("utility.jl")

# Othes
include("log.jl")

end # module
