using Base.Test
using JuMP, Ipopt, MathProgBase, Cbc, Gurobi, AmplNLWriter
using POD

# Expression Testing Instances
include("../examples/exprstest.jl")

# NLP Testing Instances
include("../examples/nlp.jl")
include("../examples/castro2m2.jl")

# Multilinear Testing Instances
include("../examples/blend029.jl")
include("../examples/multi.jl")

# Special Operator
include("../examples/div.jl")
include("../examples/circle.jl")

# Performe Tests
include("operators.jl")
include("solver.jl")
include("expression.jl")
# include("algorithm.jl")
