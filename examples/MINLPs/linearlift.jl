function basic_linear_lift(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(Alpine.Optimizer(nlp_solver=IpoptSolver(),
									mip_solver=CbcSolver(log_level=0),
									log_level=10000))
	else
		m = Model(solver)
	end

	@variable(m, x[i=1:3]>=1) # At some point if an initial value is given, keep them

    @constraint(m, (x[1]-x[2])*(3*x[2]-x[3]) >= 111)
    @NLconstraint(m, (x[1]-x[2])*(3*x[2]-x[3]) >= 111)
    @NLconstraint(m, (x[1]-x[2])^2 >= 100)
	@NLconstraint(m, (x[1]-x[2])*(3*x[2]-x[3])*(x[1]+x[3]) >= 111)
	@NLconstraint(m, (3+x[1]+x[2]+x[3])*(x[1]+x[2])^2 >= 111)

	# Next level goal
	# @NLconstraint(m, (9-x[1]-x[2]-x[3])^2 >= 111)

	@objective(m, Min, x[1]+x[2]+x[3])

	verbose && print(m)

	return m
end
