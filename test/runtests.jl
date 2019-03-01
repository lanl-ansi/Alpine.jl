using JuMP, MathOptInterface
using Ipopt, Cbc
using GLPK
using Alpine

using Test

using Random
using Pkg
alpine_dir = joinpath(dirname(pathof(Alpine)), "..")

nlp_optimizer() = Ipopt.Optimizer(print_level=0, sb="yes")
mip_optimizer() = GLPK.Optimizer()

function DefaultTestSolver(;
    nlp_optimizer=nlp_optimizer(),
    mip_optimizer=mip_optimizer(), solver_args...)
    
    solver_args_dict = Dict{Symbol,Any}()
    solver_args_dict[:nlp_optimizer] = nlp_optimizer
    solver_args_dict[:mip_optimizer] = mip_optimizer
    for v in solver_args
        solver_args_dict[v[1]] = v[2]
    end
    return solver_args_dict
end

@testset "Alpine" begin 
    include("loading.jl")
    include("convexity.jl")
end 