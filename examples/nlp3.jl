using POD, JuMP, Ipopt, Gurobi, MathProgBase

function max_cover_var_picker(m::POD.PODNonlinearModel)
	nodes = Set()
	for pair in keys(m.nonlinear_info)
		for i in pair
			@assert isa(i.args[2], Int)
			push!(nodes, i.args[2])
		end
	end
	nodes = collect(nodes)
	m.num_var_discretization_mip = length(nodes)
	m.var_discretization_mip = nodes
	return
end


function nlp3(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
								   log_level=1,
								   rel_gap=0.05,
								#    bilinear_convexhull=true,
								   presolve_bt_width_tol=1e-3,
								   presolve_bt_output_tol=1e-1,
								   presolve_bound_tightening=false,
	                               presolve_bound_tightening_algo=2,
								   discretization_var_pick_algo=0))
	else
		m = Model(solver=solver)
	end

	@variable(m, x[1:8])

	setlowerbound(x[1], 100)
	setlowerbound(x[2], 1000)
	setlowerbound(x[3], 1000)
	setlowerbound(x[4], 10)
	setlowerbound(x[5], 10)
	setlowerbound(x[6], 10)
	setlowerbound(x[7], 10)
	setlowerbound(x[8], 10)

	setupperbound(x[1], 10000)
	setupperbound(x[2], 10000)
	setupperbound(x[3], 10000)
	setupperbound(x[4], 1000)
	setupperbound(x[5], 1000)
	setupperbound(x[6], 1000)
	setupperbound(x[7], 1000)
	setupperbound(x[8], 1000)

	@constraint(m, 0.0025*(x[4]+x[6]) <= 1)
	@constraint(m, 0.0025*(x[5] - x[4] + x[7]) <= 1)
	@constraint(m, 0.01(x[8]-x[5]) <= 1)
	@NLconstraint(m, 100*x[1] - x[1]*x[6] + 833.33252*x[4] <= 83333.333)
	@NLconstraint(m, x[2]*x[4] - x[2]*x[7] - 1250*x[4] + 1250*x[5] <= 0)
	@NLconstraint(m, x[3]*x[5] - x[3]*x[8] - 2500*x[5] + 1250000 <= 0)

	@objective(m, Min, x[1]+x[2]+x[3])

	if verbose
		print(m)
	end

	return m
end
