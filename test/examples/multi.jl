# println("--------------------------------------------------------------------------")
# println("Multi4/N - exprmode 1 -> X1 * X2 * X3 * X4")
# println("Multi4/N - exprmode 2 -> (X1*X2) * (X3*X4)")
# println("Multi4/N - exprmode 3 -> (X1*X2) * X3 * X4")
# println("Multi4/N - exprmode 4 -> X1 * X2 * (X3*X4)")
# println("Multi4/N - exprmode 5 -> ((X1*X2) * X3) * X4")
# println("Multi4/N - exprmode 6 -> (X1*X2*X3) *X4")
# println("Multi4/N - exprmode 7 -> X1 * (X2 * (X3*X4))")
# println("Multi4/N - exprmode 8 -> X1 * (X2*X3) * X4")
# println("Multi4/N - exprmode 9 -> X1 * (X2*X3*X4)")
# println("Multi4/N - exprmode 10 -> X1 * ((X2*X3) * X4)")
# println("Multi4/N - exprmode 11 -> (X1 * (X2*X3)) * X4")
# println("--------------------------------------------------------------------------")
# println("Multi3/N - exprmode 1 -> X1 * X2 * X3")
# println("Multi3/N - exprmode 2 -> (X1*X2) * X3")
# println("Multi3/N - exprmode 3 -> X1 * (X2*X3)")
# println("--------------------------------------------------------------------------")
# println("                                Options")
# println("N | K | convhull | exprmode | unifrom | randomub | sos2 | presolve | delta")
# println("--------------------------------------------------------------------------")

function multi4(;verbose=false,solver=nothing, exprmode=1)

	m = Model(solver)

	Random.seed!(1)
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

function multi3(;verbose=false, solver=nothing, exprmode=1)

	m = Model(solver)

	Random.seed!(1)
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

	m = Model(solver)

	Random.seed!(1)
    @variable(m, 0.1<=x[1:2]<=rand()*10)
	@NLobjective(m, Max, x[1]*x[2])
	@constraint(m, sum(x[i] for i in 1:2) <= 2)

	if verbose
		print(m)
	end

	return m
end

function multi4N(;verbose=false, exprmode=1, solver=nothing, N=1, randomub=true)

	m = Model(solver)

	M = 1+3*N
	Random.seed!(100)
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

function multi3N(;verbose=false, randomub=nothing, solver=nothing, exprmode=1, N=1, delta=4)

	m = Model(solver)

	M = 1+2*N
	Random.seed!(100)
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

function multiKND(;verbose=false, solver=nothing, exprmode=1, randomub=true, N=1, K=2, D=1)

	m = Model(solver)

	M = K+(K-D)*(N-1)
	Random.seed!(100)
	isa(randomub, Int) ? @variable(m, 0.1<=x[1:M]<=randomub) : @variable(m, 0.1<=x[1:M]<=rand()*100)
	@NLobjective(m, Max, sum(prod(x[i+k] for k in 0:(K-1)) for i in 1:(K-D):(M-D)))
	@constraint(m, [i in 1:(K-D):(M-D)], sum(x[i+k] for k in 0:(K-1)) <= K)

	verbose && print(m)

	return m
end
