using Alpine 
using CPLEX 
using Gurobi
using JuMP 
using Ipopt
using Juniper
using Cbc
using Pavito
using Random 

include("JuMP_models.jl")
include("optimizers.jl")

# Choose underlying solvers for Alpine
mip_solver   = get_gurobi()
nlp_solver   = get_ipopt()
minlp_solver = get_juniper(mip_solver, nlp_solver)

# Global solver
const alpine = JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                            #   "minlp_solver" => minlp_solver,
                                              "nlp_solver"   => nlp_solver,  
                                              "mip_solver"   => mip_solver,
                                              "presolve_bt"  => true,
                                              "disc_ratio"   => 10)

# Try different integer values ( >=4 ) for `disc_ratio` to observe better Alpine run times, as this can be instance specific.
# Choose `presolve_bt` to `false` if you prefer the bound tightening (OBBT) presolve to be turned off. 

m = rosenbrock(solver = alpine)

JuMP.optimize!(m)
# Alpine.variable_values(m)