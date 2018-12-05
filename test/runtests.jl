using JuMP, MathProgBase
using Ipopt, Cbc, Pavito
using GLPKMathProgInterface
using POD

using Compat.Test

poddir = ""
if VERSION < v"0.7.0-"
    import Compat: round
    import Compat: undef
    import Compat: @warn
    import Compat: occursin
    import Compat: printstyled
    poddir = Pkg.dir("POD")
end

if VERSION > v"0.7.0-"
    using Random
    using Pkg
    poddir = joinpath(dirname(pathof(POD)), "..")
end

examples = readdir(joinpath(poddir, "test", "examples"))

for i in examples
    include(joinpath(poddir, "test", "examples", i))
end

pavito_solver=PavitoSolver(mip_solver=CbcSolver(logLevel=0), cont_solver=IpoptSolver(print_level=0), mip_solver_drives=false, log_level=0)

# Performe Tests
include("$(poddir)/test/solver.jl")
include("$(poddir)/test/expression.jl")
include("$(poddir)/test/algorithm.jl")
include("$(poddir)/test/utility.jl")
