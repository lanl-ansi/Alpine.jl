using JuMP, MathProgBase, Gurobi, Ipopt, POD

function circle(;verbose=false, solver=nothing, bt=false)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
                                   discretization_width_tol=1e-2,
                                   discretization_ratio=4,
								   presolve_bound_tightening=bt,
								   presolve_bound_tightening_algo=1,
								   presolve_bt_output_tol=1e-1,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

    @variable(m, 0<=x[1:2]<=2)
    @NLconstraint(m, x[1]^2 + x[2]^2 >= 2)
    @objective(m, Min, x[1]+x[2])

	if verbose
		print(m)
	end

	return m
end
