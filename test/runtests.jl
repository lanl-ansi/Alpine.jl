using JuMP, MathOptInterface
const MOI = MathOptInterface
using Ipopt, Cbc #, Pavito
using GLPK
using Alpine

using Test

using Random
using Pkg
alpine_dir = joinpath(dirname(pathof(Alpine)), "..")


#examples = readdir(joinpath(alpine_dir, "test", "examples"))

include(joinpath(alpine_dir, "test", "examples/brainpc3.jl"))
#for i in examples
#    include(joinpath(alpine_dir, "test", "examples", i))
#end

const IPOPT = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)
const IPOPT_SB = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true, "sb" => "yes")
const CBC = optimizer_with_attributes(Cbc.Optimizer, MOI.Silent() => true)

#pavito_solver=PavitoSolver(mip_solver=CbcSolver(logLevel=0), cont_solver=IpoptSolver(print_level=0, sb="yes"), mip_solver_drives=false, log_level=0)

# Perform Tests
#include("$(alpine_dir)/test/solver.jl")
include("$(alpine_dir)/test/expression.jl")
#include("$(alpine_dir)/test/algorithm.jl")
#include("$(alpine_dir)/test/utility.jl")
