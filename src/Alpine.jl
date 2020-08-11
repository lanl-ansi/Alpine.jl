__precompile__()

module Alpine

const ALPINE_DEBUG = false

using JuMP
import MathOptInterface

using LinearAlgebra
using Random
using Distributed
using Statistics
using Printf
using SparseArrays

include("const.jl")

# Engine for High-level Algorithmic Control and User-interface
include("matrix_opt_interface.jl")
include("moi_function2expr.jl")
include("solver.jl")

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
