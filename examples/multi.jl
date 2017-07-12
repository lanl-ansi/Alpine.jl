using JuMP, MathProgBase, Gurobi, Ipopt, POD

function multi(;verbose=false,solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
								   bilinear_convexhull=false,
								   presolve_bound_tightening=false,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	srand(1)
    @variable(m, 0.1<=x[1:2]<=rand()*100)
	@NLobjective(m, Max, x[1]*x[2])
	@constraint(m, sum(x[i] for i in 1:2) <= 2)

	if verbose
		print(m)
	end

	return m
end
