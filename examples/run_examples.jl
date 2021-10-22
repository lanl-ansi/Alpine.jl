using Alpine 
using CPLEX 
using Gurobi
using JuMP 
using Ipopt
using Juniper
# using Cbc
# using Pavito

include("JuMP_models.jl")
include("solver.jl")


# Choose underlying solvers for Alpine
mip_solver = get_gurobi()
nlp_solver = get_ipopt()
minlp_solver = get_juniper(mip_solver, nlp_solver)

# Global solver
const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                        #  "minlp_solver" => minlp_solver,
                                        "nlp_solver" => nlp_solver,  
                                        "mip_solver" => mip_solver,
                                        "presolve_bt" => true,
                                        "presolve_bt_max_iter" => 5,
                                        "disc_ratio" => 10)

# Try different integer values ( >=4 ) for `disc_ratio` to have better Alpine run times 
# Choose `presolve_bt` to `false` if you prefer the OBBT presolve to be turned off. 

m = nlp3(solver=alpine)

JuMP.optimize!(m)