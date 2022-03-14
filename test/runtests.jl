using JuMP, MathOptInterface
using LinearAlgebra
using Ipopt
using Cbc
using Juniper
import Pavito
using GLPK
using Alpine
using Test
using Random

alpine_dir = joinpath(dirname(pathof(Alpine)), "..")

examples = readdir(joinpath(alpine_dir, "examples", "MINLPs"))

for i in examples
    include(joinpath(alpine_dir, "examples", "MINLPs", i))
end

const MOI = MathOptInterface

const IPOPT   = optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true, "sb" => "yes", "max_iter" => 9999)
const CBC     = optimizer_with_attributes(Cbc.Optimizer, MOI.Silent() => true)
const JUNIPER = optimizer_with_attributes(Juniper.Optimizer, MOI.Silent() => true, "mip_solver" => CBC, "nl_solver" => IPOPT)
const PAVITO  = optimizer_with_attributes(Pavito.Optimizer, MOI.Silent() => true, "mip_solver" => CBC, "cont_solver" => IPOPT, "mip_solver_drives" => false)

function _build(model::JuMP.Model)
    MOI.set(model, MOI.NLPBlock(), JuMP._create_nlp_block_data(model))
    MOI.Utilities.attach_optimizer(model)
    alpine = JuMP.backend(model).optimizer.model
    Alpine.load!(alpine)
    return alpine
end

# Perform Tests
@testset "Alpine tests" begin
    include(joinpath(@__DIR__, "test_solver.jl"))
    include(joinpath(@__DIR__, "test_expression.jl"))
    include(joinpath(@__DIR__, "test_algorithm.jl"))
    include(joinpath(@__DIR__, "test_utility.jl"))
end
