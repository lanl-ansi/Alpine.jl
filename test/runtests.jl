using JuMP
using Test

import Alpine
import HiGHS
import Ipopt
import Juniper
import Pavito
import Random

const _EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples", "MINLPs")
for file in readdir(_EXAMPLES_DIR)
    include(joinpath(_EXAMPLES_DIR, file))
end

const IPOPT = MOI.OptimizerWithAttributes(
    Ipopt.Optimizer,
    MOI.Silent() => true,
    "sb" => "yes",
    "max_iter" => 9999,
)

const HIGHS = MOI.OptimizerWithAttributes(
    HiGHS.Optimizer,
    "presolve" => "on",
    "log_to_console" => false,
)

const JUNIPER = MOI.OptimizerWithAttributes(
    Juniper.Optimizer,
    MOI.Silent() => true,
    "mip_solver" => HIGHS,
    "nl_solver" => IPOPT,
)
const PAVITO = MOI.OptimizerWithAttributes(
    Pavito.Optimizer,
    MOI.Silent() => true,
    "mip_solver" => HIGHS,
    "cont_solver" => IPOPT,
    "mip_solver_drives" => false,
)

function _build(model::JuMP.Model)
    JuMP.set_optimize_hook(model, MOI.Utilities.attach_optimizer)
    JuMP.optimize!(model)
    JuMP.JuMP.set_optimize_hook(model, nothing)
    alpine = JuMP.backend(model).optimizer.model
    Alpine.load!(alpine)
    return alpine
end

# Perform Tests
@testset "Alpine tests" begin
    # include(joinpath(@__DIR__, "test_solver.jl"))
    include(joinpath(@__DIR__, "test_expression.jl"))
    # include(joinpath(@__DIR__, "test_algorithm.jl"))
    # include(joinpath(@__DIR__, "test_utility.jl"))
end
