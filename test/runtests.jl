using Base.Test
using JuMP, MathProgBase
using Ipopt, Cbc, Pajarito
using GLPKMathProgInterface
using POD

# Expression Testing Instances
include("$(Pkg.dir())/POD/test/examples/exprstest.jl")

# NLP Testing Instances
include("$(Pkg.dir())/POD/test/examples/nlp.jl")
include("$(Pkg.dir())/POD/test/examples/castrom2.jl")

# Multilinear Testing Instances
include("$(Pkg.dir())/POD/test/examples/blend029.jl")
include("$(Pkg.dir())/POD/test/examples/multi.jl")

# Special operators
include("$(Pkg.dir())/POD/test/examples/binprod.jl")
include("$(Pkg.dir())/POD/test/examples/div.jl")
include("$(Pkg.dir())/POD/test/examples/circle.jl")
include("$(Pkg.dir())/POD/test/examples/convex.jl")
include("$(Pkg.dir())/POD/test/examples/linearlift.jl")
include("$(Pkg.dir())/POD/test/examples/specialopts.jl")

# Performe Tests
# include("$(Pkg.dir())/POD/test/solver.jl")
# include("$(Pkg.dir())/POD/test/expression.jl")
include("$(Pkg.dir())/POD/test/algorithm.jl")
