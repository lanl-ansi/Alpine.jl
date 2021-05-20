using Alpine 
using CPLEX 
using JuMP 
using Ipopt

include("nlp.jl")

# MIP solver
const cplex = optimizer_with_attributes(CPLEX.Optimizer, 
                                        MOI.Silent() => true, 
                                        "CPX_PARAM_PREIND" => false) 

# Local solver
const ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                        MOI.Silent() => true, 
                                        "sb" => "yes", 
                                        "max_iter" => 9999)

# Global solver
const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                        "nlp_solver" => ipopt,
                                        # "minlp_solver" => juniper,  
                                        "mip_solver" => cplex,
                                        "presolve_bt" => true,
                                        "presolve_bt_max_iter" => 10,
                                        "disc_ratio" => 10)

# Try different integer values ( >=4 ) for `disc_ratio` to have better Alpine run times 
# Choose `presolve_bt` to `false` if you prefer the OBBT to be turned off. 

m = nlp3(solver=alpine)

JuMP.optimize!(m)