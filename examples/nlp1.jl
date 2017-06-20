using POD, JuMP, Ipopt, Gurobi, MathProgBase

function nlp1(verbose=false)

	m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0,expect_infeasible_problem="no"),
							   mip_solver=GurobiSolver(OutputFlag=0), presolve_do_bound_tightening=true, log_level=100))

    @variable(m, 1<=x[1:2]<=10)

	@NLconstraint(m, x[1]*x[2] >= 8)
	@NLobjective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5*x[1]*x[2])

	if verbose
		print(m)
	end
	return m
end
