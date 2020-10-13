function nlp1(;verbose=false,solver=nothing, convhull=false, presolve=0)

	if solver == nothing
		m = Model(Alpine.Optimizer(nlp_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
								   bilinear_convexhull=convhull,
								   monomial_convexhull=convhull,
								   presolve_bt=(presolve>0),
								   presolve_bt_algo=presolve,
								   presolve_bt_output_tol=1e-1,
								   loglevel=10000))
	else
		m = Model(solver)
	end

    @variable(m, 1<=x[1:2]<=10)

	@NLconstraint(m, x[1]*x[2] >= 8)
	@NLobjective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5*x[1]*x[2])

	return m
end

function nlp2(;verbose=false,solver=nothing, convhull=false, presolve=0)

	if solver == nothing
		m = Model(Alpine.Optimizer(nlp_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
								   bilinear_convexhull=convhull,
								   monomial_convexhull=convhull,
								   presolve_bt=(presolve>0),
								   presolve_bt_algo=presolve,
								   presolve_bt_output_tol=1e-1,
								   loglevel=10000))
	else
		m = Model(solver)
	end

	@variable(m, -500<=x[1:2]<=500)

	@NLobjective(m, Min, sum((x[i]^2 - i)^2 for i in 1:2))

	return m
end

function max_cover_var_picker(m::Alpine.Optimizer)
	nodes = Set()
	for pair in keys(m.nonconvex_terms)
		for i in pair
			@assert isa(i.args[2], Int)
			push!(nodes, i.args[2])
		end
	end
	nodes = collect(nodes)
	m.num_var_disc_mip = length(nodes)
	m.disc_vars = nodes
	return
end


function nlp3(;solver=nothing)

	m = Model(solver)

   LB = [100, 1000, 1000, 10, 10, 10, 10, 10]
   UB = [10000, 10000, 10000,  1000, 1000, 1000, 1000, 1000]

   @variable(m, LB[i] <= x[i=1:8] <= UB[i])

	@constraint(m, 0.0025*(x[4] + x[6]) <= 1)
	@constraint(m, 0.0025*(x[5] - x[4] + x[7]) <= 1)
	@constraint(m, 0.01(x[8]-x[5]) <= 1)
	@NLconstraint(m, 100*x[1] - x[1]*x[6] + 833.33252*x[4] <= 83333.333)
	@NLconstraint(m, x[2]*x[4] - x[2]*x[7] - 1250*x[4] + 1250*x[5] <= 0)
	@NLconstraint(m, x[3]*x[5] - x[3]*x[8] - 2500*x[5] + 1250000 <= 0)

	@objective(m, Min, x[1]+x[2]+x[3])

	return m
end
