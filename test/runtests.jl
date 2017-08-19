using Base.Test
using JuMP, Ipopt, MathProgBase, Cbc, Gurobi
using POD

# Expression Testing Instances
include("../examples/exprstest.jl")

# NLP Testing Instances
include("../examples/nlp.jl")
include("../examples/castrom2.jl")

# MINLP Testing Instances
include("../examples/util.jl")
include("../examples/meanvarx.jl")
include("../examples/blend.jl")

# Multilinear Testing Instances
include("../examples/multi.jl")

# Special Operator
include("../examples/div.jl")
include("../examples/circle.jl")

# Performe Tests
include("operators.jl")
include("solver.jl")
include("expression.jl")
include("algorithm.jl")
