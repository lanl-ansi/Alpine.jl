__precompile__()

module Alpine

using JuMP

import LinearAlgebra: dot, Diagonal
import Statistics
import Combinatorics

const _ALPINE_VERSION = "v0.5.2"
const ALPINE_DEBUG = false
const Alp = Alpine

include("const.jl")

# Engine for high-level algorithmic control and user-interface

include("solver_options.jl")
include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_function2expr.jl")

# Engine for expression handling
include("nlexpr.jl")
include("operators.jl")

# Main algorithm
include("main_algorithm.jl")
include("presolve.jl")
include("bounding_model.jl")
include("embedding.jl")
include("heuristics.jl")

# Convexification
include("relaxations.jl")
include("multilinear.jl")

# Model manipulation and utilities
include("variable_bounds.jl")
include("utility.jl")

# Logging
include("log.jl")

end # module
