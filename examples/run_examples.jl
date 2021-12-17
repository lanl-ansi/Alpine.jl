using Alpine 
using CPLEX 
using Gurobi
using JuMP 
using Ipopt
using Juniper

# using Cbc
# using Pavito

include("JuMP_models.jl")
include("optimizers.jl")

# Choose underlying solvers for Alpine
mip_solver = get_gurobi()
nlp_solver = get_ipopt()
minlp_solver = get_juniper(mip_solver, nlp_solver)

# Global solver
const alpine = JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                              "minlp_solver" => minlp_solver,
                                              "nlp_solver" => nlp_solver,  
                                              "mip_solver" => mip_solver,
                                              "presolve_bt" => false,
                                              "presolve_bt_max_iter" => 5,
                                              "disc_ratio" => 12)

# Try different integer values ( >=4 ) for `disc_ratio` to have better Alpine run times 
# Choose `presolve_bt` to `false` if you prefer the OBBT presolve to be turned off. 

m_alp = nlp3(solver=alpine)
JuMP.optimize!(m_alp)