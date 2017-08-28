using JuMP, MathProgBase, Gurobi, Ipopt, POD

println("--------------------------------------------------------------------------")
println("Multi4/N - exprmode 1 -> X1 * X2 * X3 * X4")
println("Multi4/N - exprmode 2 -> (X1*X2) * (X3*X4)")
println("Multi4/N - exprmode 3 -> (X1*X2) * X3 * X4")
println("Multi4/N - exprmode 4 -> X1 * X2 * (X3*X4)")
println("Multi4/N - exprmode 5 -> ((X1*X2) * X3) * X4")
println("Multi4/N - exprmode 6 -> (X1*X2*X3) *X4")
println("Multi4/N - exprmode 7 -> X1 * (X2 * (X3*X4))")
println("Multi4/N - exprmode 8 -> X1 * (X2*X3) * X4")
println("Multi4/N - exprmode 9 -> X1 * (X2*X3*X4)")
println("Multi4/N - exprmode 10 -> X1 * ((X2*X3) * X4)")
println("Multi4/N - exprmode 11 -> (X1 * (X2*X3)) * X4")
println("--------------------------------------------------------------------------")
println("Multi3/N - exprmode 1 -> X1 * X2 * X3")
println("Multi3/N - exprmode 2 -> (X1*X2) * X3")
println("Multi3/N - exprmode 3 -> X1 * (X2*X3)")
println("--------------------------------------------------------------------------")
println("                                Options")
println("N | K | convhull | exprmode | unifrom | randomub | sos2 | presolve | delta")
println("--------------------------------------------------------------------------")

function multi4(;verbose=false,solver=nothing, exprmode=1)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=0),
								   rel_gap=0.0001,
								   bilinear_convexhull=true,
								   convexhull_sweep_limit=1,
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

function multi3(;verbose=false, solver=nothing, exprmode=1, convhull=false)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=1, Presolve=0),
								   maxiter=1,
								   bilinear_convexhull=convhull,
								   discretization_add_partition_method="uniform",
								   discretization_uniform_rate=5,
								   presolve_bound_tightening=false,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	srand(1)
	ub = rand()
    @variable(m, 0.1<=x[1:3]<=ub*1000)
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
								   mip_solver=GurobiSolver(OutputFlag=0),
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

function multi4N(;verbose=false, solver=nothing, N=1, convhull=false, exprmode=1, uniform=-1, randomub=true, sos2=true, delta=4, presolve=0)

	if solver == nothing
		if uniform >0
			m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									   mip_solver=GurobiSolver(OutputFlag=1),
									   bilinear_convexhull=convhull,
									   convexhull_use_sos2=sos2,
									   discretization_add_partition_method="uniform",
									   discretization_uniform_rate=uniform,
									   maxiter=1,
									   log_level=100))
		else
			m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									   mip_solver=GurobiSolver(OutputFlag=0),
									   discretization_ratio=delta,
									   convexhull_use_sos2=sos2,
									   bilinear_convexhull=convhull,
									   presolve_bound_tightening=(presolve>0),
									   presolve_bound_tightening_algo=presolve,
									   log_level=100))
		end
	else
		m = Model(solver=solver)
	end

	M = 1+3*N
	srand(100)
	isa(randomub, Bool) ? @variable(m, 0.1 <= x[1:M] <= 100*rand()) : @variable(m, 0.1 <= x[1:M] <= randomub)

	if exprmode == 1
		@NLobjective(m, Max, sum(x[i]*x[i+1]*x[i+2]*x[i+3] for i in 1:3:(M-1)))
	elseif exprmode == 2
		@NLobjective(m, Max, sum((x[i]*x[i+1])*(x[i+2]*x[i+3]) for i in 1:3:(M-1)))
	elseif exprmode == 3
		@NLobjective(m, Max, sum((x[i]*x[i+1])*x[i+2]*x[i+3] for i in 1:3:(M-1)))
	elseif exprmode == 4
		@NLobjective(m, Max, sum(x[i]*x[i+1]*(x[i+2]*x[i+3]) for i in 1:3:(M-1)))
	elseif exprmode == 5
		@NLobjective(m, Max, sum((((x[i]*x[i+1])*x[i+2])*x[i+3]) for i in 1:3:(M-1)))
	elseif exprmode == 6
		@NLobjective(m, Max, sum((x[i]*x[i+1]*x[i+2])*x[i+3] for i in 1:3:(M-1)))
	elseif exprmode == 7
		@NLobjective(m, Max, sum(x[i]*(x[i+1]*(x[i+2]*x[i+3])) for i in 1:3:(M-1)))
	elseif exprmode == 8
		@NLobjective(m, Max, sum(x[i]*(x[i+1]*x[i+2])*x[i+3] for i in 1:3:(M-1)))
	elseif exprmode == 9
		@NLobjective(m, Max, sum(x[i]*(x[i+1]*x[i+2]*x[i+3]) for i in 1:3:(M-1)))
	elseif exprmode == 10
		@NLobjective(m, Max, sum(x[i]*((x[i+1]*x[i+2])*x[i+3]) for i in 1:3:(M-1)))
	elseif exprmode == 11
		@NLobjective(m, Max, sum((x[i]*(x[i+1]*x[i+2]))*x[i+3] for i in 1:3:(M-1)))
	else
		error("exprmode argument only taks from 1-7")
	end

	@constraint(m, [i in 1:3:(M-1)], x[i]+x[i+1]+x[i+2]+x[i+3]<=4)

	if verbose
		print(m)
	end

	return m
end

function multi3N(;verbose=false, solver=nothing, exprmode=1, convhull=false, uniform=-1, randomub=true, N=1, delta=4, sos2=false, presolve=0)

	if solver == nothing
		if uniform > 0.0
			m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									   mip_solver=GurobiSolver(OutputFlag=1),
									   convexhull_use_sos2=sos2,
									   maxiter=1,
									   bilinear_convexhull=convhull,
									   discretization_add_partition_method="uniform",
									   discretization_uniform_rate=uniform,
									   presolve_bound_tightening=false,
									   log_level=1))
		else
			m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									   mip_solver=GurobiSolver(OutputFlag=0),
									   convexhull_use_sos2=sos2,
									   discretization_ratio=delta,
									   bilinear_convexhull=convhull,
									   presolve_bound_tightening=(presolve>0),
									   presolve_bound_tightening_algo=presolve,
									   discretization_var_pick_algo=0,
									   log_level=1))
		end
	else
		m = Model(solver=solver)
	end

	M = 1+2*N
	srand(100)
    isa(randomub, Int) ? @variable(m, 0.1<=x[1:M]<=randomub) : @variable(m, 0.1<=x[1:M]<=rand()*100)
	if exprmode == 1
		@NLobjective(m, Max, sum(x[i]*x[i+1]*x[i+2] for i in 1:2:(M-1)))
	elseif exprmode == 2
		@NLobjective(m, Max, sum((x[i]*x[i+1])*x[i+2] for i in 1:2:(M-1)))
	elseif exprmode == 3
		@NLobjective(m, Max, sum(x[i]*(x[i+1]*x[i+2]) for i in 1:2:(M-1)))
	else
		error("exprmode argument only taks from 1-7")
	end

	@constraint(m, [i in 1:2:(M-1)], x[i] + x[i+1] + x[i+2] <= 3)

	if verbose
		print(m)
	end

	return m
end

function multiKND(;verbose=false, solver=nothing, exprmode=1, convhull=false, uniform=-1, sos2=false, facet=false, randomub=true, N=1, K=2, D=1, delta=4, presolve=0)

	if solver == nothing
		if uniform > 0.0
			m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									   mip_solver=GurobiSolver(OutputFlag=1),
									   maxiter=1,
									   bilinear_convexhull=convhull,
									   convexhull_use_sos2=sos2,
									   convexhull_use_facet=facet,
									   discretization_add_partition_method="uniform",
									   discretization_uniform_rate=uniform,
									   presolve_bound_tightening=false,
									   log_level=100))
		else
			m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									   mip_solver=GurobiSolver(OutputFlag=1),
									   convexhull_use_sos2=sos2,
									   convexhull_use_facet=facet,
									   discretization_ratio=delta,
									   bilinear_convexhull=convhull,
									   presolve_bound_tightening=(presolve>0),
									   presolve_bound_tightening=presolve,
									   log_level=100))
		end
	else
		m = Model(solver=solver)
	end

	M = K+(K-D)*(N-1)
	srand(100)
	isa(randomub, Int) ? @variable(m, 0.1<=x[1:M]<=randomub) : @variable(m, 0.1<=x[1:M]<=rand()*100)
	@NLobjective(m, Max, sum(prod(x[i+k] for k in 0:(K-1)) for i in 1:(K-D):(M-D)))
	@constraint(m, [i in 1:(K-D):(M-D)], sum(x[i+k] for k in 0:(K-1)) <= K)

	verbose && print(m)

	return m
end
