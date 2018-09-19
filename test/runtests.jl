using Base.Test
using JuMP, MathProgBase
using Ipopt, Cbc, Pavito
using GLPKMathProgInterface
using POD

Pkg.pin("Cbc", v"0.3.2")

poddir = Pkg.dir("POD")

examples = readdir(joinpath(Pkg.dir("POD"), "test", "examples"))
for i in examples
    include(joinpath(Pkg.dir("POD"), "test", "examples", i))
end

pavito_solver=PavitoSolver(mip_solver=CbcSolver(logLevel=0), cont_solver=IpoptSolver(print_level=0), mip_solver_drives=false, log_level=0)

# Performe Tests
include("$(poddir)/test/solver.jl")
include("$(poddir)/test/expression.jl")
include("$(poddir)/test/algorithm.jl")
include("$(poddir)/test/utility.jl")
