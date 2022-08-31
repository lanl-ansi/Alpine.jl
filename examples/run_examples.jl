# - Necessary - #
using Alpine
using JuMP
using Gurobi
using Ipopt
using Juniper

# - Additional - #
# using CPLEX
# using HiGHS
# using Pavito

include("JuMP_models.jl")
include("optimizers.jl")

# Choose underlying solvers for Alpine
nlp_solver = get_ipopt() #local solver
mip_solver = get_gurobi()
minlp_solver = get_juniper(mip_solver, nlp_solver) #local solver

#= Global solver
 Hints: 
 => Try different integer values (>=4) for `disc_ratio` to potentially observe 
    better Alpine run times (can be instance specific).
 => Choose `presolve_bt` to `false` if you prefer the bound tightening (OBBT) presolve to be turned off. 
 => If you prefer to use Alpine for only OBBT presolve, without any paritioning applied to the 
    nonlinear terms, include option "apply_partitioning" below and set it to false. 
=#

const alpine = JuMP.optimizer_with_attributes(
    Alpine.Optimizer,
    # "minlp_solver" => minlp_solver,
    "nlp_solver" => nlp_solver,
    "mip_solver" => mip_solver,
    "presolve_bt" => true,
    "disc_ratio" => 10,
)

m = nlp3(solver = alpine)

JuMP.optimize!(m)
# Alpine.variable_values(m)
