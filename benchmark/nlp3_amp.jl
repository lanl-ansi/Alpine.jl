using JuMP, MathProgBase
using Gurobi, Ipopt
using POD

include("../examples/nlp3.jl")
benchmark_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
		  					mip_solver=GurobiSolver(OutputFlag=0),
						   	maxiter=3, log_level=1)

m = nlp3(solver=benchmark_solver)
status = solve(m)
println("Final POD Status is $status")
