using JuMP, MathProgBase
using Gurobi, Ipopt
using POD

include("../examples/nlp1.jl")
benchmark_solver = PODSolver(nlp_local_solver=IpoptSolver(),
							   mip_solver=GurobiSolver(OutputFlag=0),
							   log_level=1,
							   presolve_bt_width_tol=1e-3,
							   presolve_bt_output_tol=1e-1,
							   presolve_perform_bound_tightening=true,
							   presolve_bound_tightening_algo=1)

m = nlp1(solver=benchmark_solver)
status = solve(m)
println("Final POD Status is $status")
