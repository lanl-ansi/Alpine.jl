using JuMP, MathProgBase
using Ipopt, Cbc, Pavito
using GLPKMathProgInterface
using Alpine

using Compat.Test

alpine_dir = ""
if VERSION < v"0.7.0-"
    import Compat: round
    import Compat: undef
    import Compat: @warn
    import Compat: occursin
    import Compat: printstyled
    alpine_dir = Pkg.dir("Alpine")
end

if VERSION > v"0.7.0-"
    using Random
    using Pkg
    alpine_dir = joinpath(dirname(pathof(Alpine)), "..")
end

examples = readdir(joinpath(alpine_dir, "test", "examples"))

for i in examples
    include(joinpath(alpine_dir, "test", "examples", i))
end

pavito_solver=PavitoSolver(mip_solver=CbcSolver(logLevel=0), cont_solver=IpoptSolver(print_level=0), mip_solver_drives=false, log_level=0)
