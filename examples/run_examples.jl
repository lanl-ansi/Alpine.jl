using Alpine 
using CPLEX 
using Gurobi
using JuMP 
using Ipopt
using Juniper
using Cbc
using Pavito

include("JuMP_models.jl")
include("optimizers.jl")

# Choose underlying solvers for Alpine
cbc = get_cbc()
ipopt = get_ipopt()
gurobi = get_gurobi()
mip_solver = get_pavito(gurobi, ipopt)
nlp_solver = get_ipopt()
minlp_solver = get_juniper(mip_solver, nlp_solver)

# Global solver
const alpine = JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                            #  "minlp_solver" => minlp_solver,
                                              "nlp_solver" => nlp_solver,  
                                              "mip_solver" => mip_solver,
                                              "presolve_bt" => false,
                                              "presolve_bt_max_iter" => 5,
                                              "disc_ratio" => 8,
                                              "disc_abs_width_tol" => 1e-2,
                                              "max_iter" => 20)

# Try different integer values ( >=4 ) for `disc_ratio` to have better Alpine run times 
# Choose `presolve_bt` to `false` if you prefer the OBBT presolve to be turned off. 

m = circle(solver=alpine)

JuMP.optimize!(m)