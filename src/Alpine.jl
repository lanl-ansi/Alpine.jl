__precompile__()

module Alpine

using JuMP

import LinearAlgebra: dot, Diagonal
import Statistics

const ALPINE_DEBUG = false
const Alp = Alpine

include("const.jl")

# Engine for high-level algorithmic control and user-interface
include("matrix_opt_interface.jl")
include("moi_function2expr.jl")
include("solver.jl")

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
