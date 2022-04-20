__precompile__()

module Alpine

using JuMP

import LinearAlgebra: dot, Diagonal
import Statistics

const ALPINE_DEBUG = false
const Alp = Alpine

include("const.jl")

# Engine for high-level algorithmic control and user-interface

include("solver.jl")
include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_function2expr.jl")

# Engine for expression handling
include("nlexpr.jl")
include("operators.jl")

# Main algorithm
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

# Logging
include("log.jl")

end # module
