# Contains a basic model with various expressions for testing
using POD, JuMP, Ipopt, Cbc, MathProgBase

function basic_linear_lift(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									mip_solver=CbcSolver(),
									log_level=0))
	else
		m = Model(solver=solver)
	end

	@variable(m, x[i=1:6]>=1) # At some point if an initial value is given, keep them

    @constraint(m, (x[1]-x[2])*(3*x[2]-x[3]) >= 111)
    @NLconstraint(m, (x[1]-x[2])*(3*x[2]-x[3]) >= 111)
    @NLconstraint(m, (x[1]-x[2])^2 >= 100)
    @NLconstraint(m, (x[1]-x[2])*(3*x[2]-x[3])*(x[4]+x[6]) >= 111)

	@objective(m, Min, x[1]+x[2]+x[3])

	if verbose
		print(m)
	end

	return m
end
