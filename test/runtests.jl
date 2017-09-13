using Base.Test
using JuMP, Ipopt, MathProgBase, Cbc, Gurobi, AmplNLWriter
using POD

# Expression Testing Instances
include("$(Pkg.dir())/POD/examples/exprstest.jl")

# NLP Testing Instances
include("$(Pkg.dir())/POD/examples/nlp.jl")
include("$(Pkg.dir())/POD/examples/castro2m2.jl")

# Multilinear Testing Instances
include("$(Pkg.dir())/POD/examples/blend029.jl")
include("$(Pkg.dir())/POD/examples/multi.jl")

# Special Operator
include("$(Pkg.dir())/POD/examples/div.jl")
include("$(Pkg.dir())/POD/examples/circle.jl")
include("$(Pkg.dir())/POD/examples/convex.jl")
include("$(Pkg.dir())/POD/examples/linearlift.jl")

# Performe Tests
include("$(Pkg.dir())/POD/test/solver.jl")
include("$(Pkg.dir())/POD/test/expression.jl")
include("$(Pkg.dir())/POD/test/algorithm.jl")
