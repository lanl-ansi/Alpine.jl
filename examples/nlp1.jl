using POD, JuMP, Ipopt, MathProgBase 

function pod_example_nlp1(verbose=false)

	m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0)))
	@variable(m, 1<=x[1:2]<=10)

	@NLconstraint(m, x[1]*x[2] >= 8)
	@NLobjective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5*x[1]*x[2])

	if verbose
		print(m)
	end
	return m
end
