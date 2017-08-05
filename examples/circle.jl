using JuMP, MathProgBase, Gurobi, Ipopt, POD

function circle(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
                                   discretization_width_tol=1e-2,
                                   discretization_ratio=4,
								   presolve_bound_tightening=false,
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

function circleN(;verbose=false, solver=nothing, convhull=false, N=2)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=1),
								   monomial_convexhull=convhull,
								   discretization_width_tol=1e-2,
								   maxiter=1,
								   discretization_add_partition_method="uniform",
								   discretization_uniform_rate=100,
								   presolve_bound_tightening=false,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

    @variable(m, 0<=x[1:N]<=N)
    @NLconstraint(m, sum(x[i]^2 for i in 1:N) >= N)
    @objective(m, Min, sum(x))

	if verbose
		print(m)
	end

	return m
end
