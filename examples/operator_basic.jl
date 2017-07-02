using JuMP, MathProgBase, Ipopt, Cbc, JLD, POD

function operator_basic(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(),
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	@variable(m, x[1:4]>=0)

	@NLconstraint(m, x[1]^2 >= 1)
	@NLconstraint(m, x[1]*x[2] <= 1)
	@NLconstraint(m, x[1]^2 + x[2]*x[3] <= 1)

	@NLconstraint(m, x[1] * (x[2]*x[3]) >= 1)
	@NLconstraint(m, x[1]^2 * (x[2]^2 * x[3]^2) <= 1)

	@NLconstraint(m, (x[1] * x[2]) * x[3] >= 1)
	@NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

	@NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
	@NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2 <= 1)
	@NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
	@NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

	@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 1)
	@NLconstraint(m, ((x[1]^2*x[2])*x[3])*x[4] <= 1)
	@NLconstraint(m, ((x[1]*x[2]^2)*x[3])*x[4] >= 1)
	@NLconstraint(m, ((x[1]*x[2])*x[3]^2)*x[4] <= 1)
	@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, ((x[1]^2*x[2]^2)*x[3]^2)*x[4]^2 <=1)

	@NLconstraint(m, x[1]*(x[2]*(x[3]*x[4])) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*(x[3]*x[4])) <= 1)
	@NLconstraint(m, x[1]*(x[2]^2*(x[3]*x[4])) >= 1)
	@NLconstraint(m, x[1]*(x[2]*(x[3]^2*x[4])) <= 1)
	@NLconstraint(m, x[1]*(x[2]*(x[3]*x[4]^2)) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]^2*(x[3]^2*x[4]^2)) <= 1)

	@NLconstraint(m, x[1]*x[2]*x[3] >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3] >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3] >= 1)
	@NLconstraint(m, x[1]*x[2]*x[3]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*x[3] >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2 >= 1)

	@NLconstraint(m, (x[1]*x[2])*(x[3]*x[4]) >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]) >= 1)
	@NLconstraint(m, (x[1]*x[2]^2)*(x[3]*x[4]) >= 1)
	@NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]) >= 1)
	@NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]^2) >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]^2) >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*(x[3]^2*x[4]) >= 1)

	@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3]*x[4] >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3]*x[4] >= 1)
	@NLconstraint(m, x[1]*x[2]*x[3]^2*x[4]^2 >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3]^2*x[4] >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3]*x[4]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4]^2 >= 1)

	@NLconstraint(m, (x[1]*x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, (x[1]^2*x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]^2*x[3])*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]*x[3]^2)*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, (x[1]^2*x[2]^2*x[3]^2)*x[4]^2 >= 1)

	@NLconstraint(m, x[1]*(x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]^2*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, x[1]*(x[2]^2*x[3]^2)*x[4] >= 1)
	@NLconstraint(m, x[1]^2*(x[2]^2*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4]^2 >= 1)

	@NLconstraint(m, x[1]*(x[2]*x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*(x[2]^2*x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]^2*x[4]) >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]*x[4]^2) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3]^2*x[4]) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]^2*x[3]^2*x[4]^2) >= 1)

	@NLconstraint(m, x[1]*x[2]*(x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]^2*x[2]*(x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*x[2]^2*(x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*x[2]*(x[3]^2*x[4]) >= 1)
	@NLconstraint(m, x[1]*x[2]*(x[3]*x[4]^2) >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*(x[3]^2*x[4]^2) >= 1)

	@NLconstraint(m, (x[1]*x[2])*x[3]*x[4] >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*x[3]*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]^2)*x[3]*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2])*x[3]*x[4]^2 >= 1)
	@NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4]^2 >= 1)
	@NLconstraint(m, (x[1]^2*x[2]^2)*x[3]^2*x[4]^2 >= 1)

	if verbose
		print(m)
	end

	return m
end
