using Alpine 
using CPLEX 
using JuMP 
using Ipopt
using Cbc
using Pavito

include("nlp.jl")

# MIP solver
const cplex = optimizer_with_attributes(CPLEX.Optimizer, 
                                        MOI.Silent() => true, 
                                        "CPX_PARAM_PREIND" => false) 

const cbc = optimizer_with_attributes(Cbc.Optimizer, 
                                      MOI.Silent() => true) 
                                     
# Local solver
const ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                        MOI.Silent() => true, 
                                        "sb" => "yes", 
                                        "max_iter" => 9999)

# Convex MINLP solver
const pavito = optimizer_with_attributes(Pavito.Optimizer, 
                                        MOI.Silent() => true,
                                        "mip_solver" => cbc, 
                                        "cont_solver" => ipopt, 
                                        "mip_solver_drives" => false)

# Global solver
const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                        # "minlp_solver" => juniper,
                                        "nlp_solver" => ipopt,  
                                        "mip_solver" => cplex,
                                        "presolve_bt" => true,
                                        "presolve_bt_max_iter" => 5,
                                        "disc_ratio" => 12)

# Try different integer values ( >=4 ) for `disc_ratio` to have better Alpine run times 
# Choose `presolve_bt` to `false` if you prefer the OBBT to be turned off. 

m = nlp3(solver=alpine)

JuMP.optimize!(m)