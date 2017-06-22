benchmark_start = time()
include("../examples/nlp1.jl")
m = nlp1()
PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
						   mip_solver=GurobiSolver(OutputFlag=0),
						   log_level=100, maxiter=4,
						   presolve_bt_width_tolerance=1e-3,
						   presolve_bt_output_tolerance=1e-1,
						   presolve_perform_bound_tightening=true,
						   presolve_bound_tightening_algo=2,
						   discretization_var_pick_algo=max_cover_var_picker)
@time status = solve(m)
benchmark_time = time() - benchmark_start
println("Algorithm Final Status: $status")
println("Total time: $(benchmark_time)s")
