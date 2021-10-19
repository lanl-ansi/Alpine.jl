using Alpine 
using CPLEX 
using Gurobi
using JuMP 
using Ipopt
using Cbc
using Pavito
using Juniper

include("JuMP_models.jl")

# MIP solver - commercial 
const gurobi = optimizer_with_attributes(Gurobi.Optimizer, 
                                        MOI.Silent() => true, 
                                        "Presolve" => 1) 

const cplex = optimizer_with_attributes(CPLEX.Optimizer, 
                                        MOI.Silent() => true, 
                                        "CPX_PARAM_PREIND" => false) 

# MIP solver - open-source
const cbc = optimizer_with_attributes(Cbc.Optimizer, 
                                      MOI.Silent() => true) 

# Choose the MIP solver here (for Alpine/Juniper)
mip_solver = gurobi  

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

const juniper = optimizer_with_attributes(Juniper.Optimizer, 
                                          MOI.Silent() => true, 
                                          "mip_solver" => mip_solver, 
                                          "nl_solver" => ipopt)

# Global solver
const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                         "minlp_solver" => juniper,
                                         "nlp_solver" => ipopt,  
                                         "mip_solver" => mip_solver,
                                         "presolve_bt" => true,
                                         "presolve_bt_max_iter" => 5,
                                         "disc_ratio" => 10)

# Try different integer values ( >=4 ) for `disc_ratio` to have better Alpine run times 
# Choose `presolve_bt` to `false` if you prefer the OBBT presolve to be turned off. 

m = nlp3(solver=alpine)

JuMP.optimize!(m)