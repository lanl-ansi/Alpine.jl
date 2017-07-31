using JuMP, MathProgBase, Gurobi, Ipopt, POD

function div(;verbose=false,solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
								   presolve_bound_tightening=true,
								   presolve_bound_tightening_algo=2,
								   presolve_bt_output_tol=1e-1,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

    @variable(m, 1<=x[1:2]<=10)

	@NLobjective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5/5*x[1]*x[2])

    @constraint(m, x[1]*(-1*3/5*5) >= 0)
    @constraint(m, x[2]*5/15 >= 0)
    @constraint(m, x[2]*(5*15/15) >= 0)
    @constraint(m, x[2]*(123.23*1231/-1231.1231) >= 0)
	@NLconstraint(m, (3/5)*x[1]*(60/60)*x[2] >= 8)

	if verbose
		print(m)
	end
	return m
end
