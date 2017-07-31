using Base.Test
using JuMP, Ipopt, MathProgBase, Cbc, Gurobi
using POD

# Expression Testing Instances
include("../examples/operator_c.jl")
include("../examples/exprstest.jl")

# NLP Testing Instances
include("../examples/nlp1.jl")
include("../examples/nlp2.jl")
include("../examples/nlp3.jl")
include("../examples/castro02m2.jl")
include("../examples/castro03m2.jl")

# MINLP Testing Instances
include("../examples/util.jl")
include("../examples/meanvarx.jl")
include("../examples/blend029.jl")

# Multilinear Testing Instances
include("../examples/multi.jl")

# Special Operator
include("../examples/div.jl")

# Performe Tests
include("operators.jl")
include("solver.jl")
include("expression.jl")
include("algorithm.jl")
