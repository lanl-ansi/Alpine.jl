function pod_example_bi2(verbose=false)

	m = Model(solver=IpoptSolver())
	@variable(m, 1<=x[i=1:10]<=5)

	@NLconstraint(m, sum((i+100)*x[i]^2 for i=1:10) >= 128)
	@NLconstraint(m, x[1]*2*x[2]*(89-13) - 4*(1-5)*x[1]*x[5] >= 222)
	@constraint(m, sum(i*x[i] for i=1:10) <= 999)
	@NLconstraint(m, sum((i-5)*x[i]^2 for i=1:10) <= 256)
	partA = randperm(10)
	partB = randperm(10)
	for i in 1:10
		@NLconstraint(m, x[partA[i]]*x[partB[i]] >= 5)
	end

	@NLobjective(m, Min, sum(x[i]*x[i+1] for i=1:9))
    if verbose
        print(m)
    end
	return m
end
