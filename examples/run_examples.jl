using Alpine 
using CPLEX 
using Gurobi
using JuMP 
using Ipopt
using Juniper
using Cbc
using Pavito
using BARON
import Random

include("JuMP_models.jl")
include("optimizers.jl")

# Choose underlying solvers for Alpine
cbc = get_cbc()
ipopt = get_ipopt()
gurobi = get_gurobi()
juniper = get_juniper(gurobi, ipopt)

const IPOPT   = MOI.OptimizerWithAttributes(Ipopt.Optimizer, MOI.Silent() => true, "sb" => "yes", "max_iter" => 9999)
const CBC     = MOI.OptimizerWithAttributes(Cbc.Optimizer, MOI.Silent() => true)
const PAVITO  = MOI.OptimizerWithAttributes(Pavito.Optimizer, MOI.Silent() => true, "mip_solver" => CBC, "cont_solver" => IPOPT, "mip_solver_drives" => false)
const JUNIPER = MOI.OptimizerWithAttributes(Juniper.Optimizer, MOI.Silent() => true, "mip_solver" => CBC, "nl_solver" => IPOPT)

const baron = optimizer_with_attributes(BARON.Optimizer, 
                                       "maxtime" => 3600, 
                                       "CplexLibName" => "/Users/harsha/Documents/CPLEX_20_1/cplex/bin/x86-64_osx/libcplex2010.jnilib")


mip_solver = get_gurobi()
nlp_solver = get_ipopt()
# minlp_solver = get_juniper(mip_solver, nlp_solver)

# Global solver
const alpine = JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                            #  "minlp_solver" => minlp_solver,
                                              "nlp_solver" => nlp_solver,  
                                              "mip_solver" => mip_solver,
                                              "presolve_bt" => true,
                                              "presolve_bt_max_iter" => 5,
                                              "disc_ratio" => 5)

# Try different integer values ( >=4 ) for `disc_ratio` to observe better Alpine run times, as this can be instance specific.
# Choose `presolve_bt` to `false` if you prefer the bound tightening (OBBT) presolve to be turned off. 

m = nlp3(solver = alpine)

JuMP.optimize!(m)
# Alpine.variable_values(m)