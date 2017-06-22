using POD, JuMP, Ipopt, Gurobi, MathProgBase

function nlp1(verbose=false)

	m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
							   mip_solver=GurobiSolver(OutputFlag=0),
							   presolve_perform_bound_tightening=false,
							   presolve_bound_tightening_algo=1,
							   presolve_bt_output_tolerance=1e-1,
							   log_level=1))

    @variable(m, 1<=x[1:2]<=10)

	@NLconstraint(m, x[1]*x[2] >= 8)
	@NLobjective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5*x[1]*x[2])

	if verbose
		print(m)
	end
	return m
end
