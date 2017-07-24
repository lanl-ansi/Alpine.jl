using JuMP, MathProgBase, CPLEX, Ipopt, POD

println("--------------------------------------------")
println("Multi4 - exprmode 1 -> X1 * X2 * X3 * X4")
println("Multi4 - exprmode 2 -> (X1*X2) * (X3*X4)")
println("Multi4 - exprmode 3 -> (X1*X2) * X3 * X4")
println("Multi4 - exprmode 4 -> X1 * X2 * (X3*X4)")
println("Multi4 - exprmode 5 -> ((X1*X2) * X3) * X4")
println("Multi4 - exprmode 6 -> (X1*X2*X3) *X4")
println("Multi4 - exprmode 7 -> X1 * (X2 * (X3*X4))")
println("Multi4 - exprmode 8 -> X1 * (X2*X3) * X4")
println("Multi4 - exprmode 9 -> X1 * (X2*X3*X4)")
println("Multi4 - exprmode 10 -> X1 * ((X2*X3) * X4)")
println("Multi4 - exprmode 11 -> (X1 * (X2*X3)) * X4")
println("--------------------------------------------")

function pick_my_var(m::Any)
	exclude = 3		# Give the variabel you wish to exclude from discretization set
	# Perform a max-cover locally and exclude the target term
	nodes = Set()
	for pair in keys(m.nonlinear_info)
		# Assumption Max cover is always safe
		for i in pair
			@assert isa(i.args[2], Int)
			push!(nodes, i.args[2])
		end
	end
	delete!(nodes, exclude)
	nodes = collect(nodes)
	m.num_var_discretization_mip = length(nodes)
	m.var_discretization_mip = nodes

	return
end

function multi4(;verbose=false,solver=nothing, exprmode=1)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=CplexSolver(CPX_PARAM_SCRIND=0),
								   rel_gap=0.0001,
								   bilinear_convexhull=true,
								   convhull_sweep_limit=1,
								#    discretization_var_pick_algo=pick_my_var,
								   presolve_bound_tightening=false,
								   presolve_track_time=true,
								   presolve_bound_tightening_algo=1,
								   log_level=100))
	else
		m = Model(solver=solver)
	end

	srand(1)
    @variable(m, 0.1<=x[1:4]<=rand()*100)
	if exprmode == 1
		@NLobjective(m, Max, x[1]*x[2]*x[3]*x[4])
	elseif exprmode == 2
		@NLobjective(m, Max, (x[1]*x[2])*(x[3]*x[4]))
	elseif exprmode == 3
		@NLobjective(m, Max, (x[1]*x[2])*x[3]*x[4])
	elseif exprmode == 4
		@NLobjective(m, Max, x[1]*x[2]*(x[3]*x[4]))
	elseif exprmode == 5
		@NLobjective(m, Max, (((x[1]*x[2])*x[3])*x[4]))
	elseif exprmode == 6
		@NLobjective(m, Max, (x[1]*x[2]*x[3])*x[4])
	elseif exprmode == 7
    	@NLobjective(m, Max, x[1]*(x[2]*(x[3]*x[4])))
	elseif exprmode == 8
		@NLobjective(m, Max, x[1]*(x[2]*x[3])*x[4])
	elseif exprmode == 9
		@NLobjective(m, Max, x[1]*(x[2]*x[3]*x[4]))
	elseif exprmode == 10
		@NLobjective(m, Max, x[1]*((x[2]*x[3])*x[4]))
	elseif exprmode == 11
		@NLobjective(m, Max, (x[1]*(x[2]*x[3]))*x[4])
	else
		error("exprmode argument only taks from 1-7")
	end
	@constraint(m, sum(x[i] for i in 1:4) <= 4)

	if verbose
		print(m)
	end

	return m
end

println("--------------------------------------------")
println("Multi3 - exprmode 1 -> X1 * X2 * X3")
println("Multi3 - exprmode 2 -> (X1*X2) * X3")
println("Multi3 - exprmode 3 -> X1 * (X2*X3)")
println("--------------------------------------------")

function multi3(;verbose=false,solver=nothing, exprmode=1)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=CplexSolver(CPX_PARAM_SCRIND=0),
								   rel_gap=0.001,
								   bilinear_convexhull=true,
								#    discretization_add_partition_method="uniform",
								#    discretization_uniform_rate=3,
								   presolve_bound_tightening=false,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	srand(1)
	ub = rand()
    @variable(m, 0.1<=x[1:3]<=ub*100)
	if exprmode == 1
		@NLobjective(m, Max, x[1]*x[2]*x[3])
	elseif exprmode == 2
		@NLobjective(m, Max, (x[1]*x[2])*x[3])
	elseif exprmode == 3
		@NLobjective(m, Max, x[1]*(x[2]*x[3]))
	else
		error("exprmode argument only taks from 1-7")
	end

	@constraint(m, sum(x[i] for i in 1:3) <= 3)

	if verbose
		print(m)
	end

	return m
end


function multi2(;verbose=false,solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=CplexSolver(CPX_PARAM_SCRIND=0),
								   bilinear_convexhull=true,
								   discretization_add_partition_method="uniform",
								   presolve_bound_tightening=false,
								   log_level=100))
	else
		m = Model(solver=solver)
	end

	srand(1)
    @variable(m, 0.1<=x[1:2]<=rand()*10)
	@NLobjective(m, Max, x[1]*x[2])
	@constraint(m, sum(x[i] for i in 1:2) <= 2)

	if verbose
		print(m)
	end

	return m
end
