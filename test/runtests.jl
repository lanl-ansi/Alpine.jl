using Base.Test
using JuMP, MathProgBase
using Ipopt, Cbc, Pajarito
using GLPKMathProgInterface
using POD

poddir = Pkg.dir("POD")

# Expression Testing Instances
include("$(poddir)/test/examples/exprstest.jl")

# NLP Testing Instances
include("$(poddir)/test/examples/nlp.jl")
include("$(poddir)/test/examples/castrom2.jl")

# Multilinear Testing Instances
include("$(poddir)/test/examples/blend029.jl")
include("$(poddir)/test/examples/multi.jl")

# Special operators
include("$(poddir)/test/examples/binprod.jl")
include("$(poddir)/test/examples/div.jl")
include("$(poddir)/test/examples/circle.jl")
include("$(poddir)/test/examples/convex.jl")
include("$(poddir)/test/examples/linearlift.jl")
include("$(poddir)/test/examples/specialopts.jl")

# Performe Tests
include("$(poddir)/test/solver.jl")
include("$(poddir)/test/expression.jl")
include("$(poddir)/test/algorithm.jl")
include("$(poddir)/test/utility.jl")
