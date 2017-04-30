# Example Model nlp2
function pod_example_nlp2(verbose=false)

	info("This model is currently not suitable for expression operations...")

	m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0)))
	@variable(m, -500<=x[1:2]<=500)

	@NLobjective(m, Min, sum((x[1]^2 - i)^2 for i in 1:2))

	if verbose
		print(m)
	end

	return m
end
