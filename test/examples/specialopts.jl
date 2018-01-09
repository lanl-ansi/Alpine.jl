function specialopts(;verbose=false, solver=nothing)

	m = Model(solver=solver)

	@variable(m, x[i=1:6]) # At some point if an initial value is given, keep them

    @NLconstraint(m, sin(x[1]) + cos(x[2]) >= 1)
    @NLconstraint(m, sin(x[1]) * x[2] >= 1)
    @NLconstraint(m, sin(x[1]) * cos(x[2]) >= 1)
    #
    @NLconstraint(m, sin(x[1]*x[2]) >= 0.5)
    @NLconstraint(m, cos(x[1]+x[2]) >= 0.5)
    @NLconstraint(m, cos(x[1]/2) >= 0.5)
    @NLconstraint(m, cos(x[1]*x[2]/5) >= 0.5)
    @NLconstraint(m, sin(x[1]/5+x[2]/5) >= 0.5)
    #

	return m
end

function bpml(;verbose=false, solver=nothing)

	m = Model(solver=solver)

	@variable(m, x[1:5], Bin)
	@variable(m, y[1:5]>=0)

	@NLconstraint(m, x[1]*y[1] >= 10)
	@NLconstraint(m, x[1]*x[2]*y[1] >= 10)
	@NLconstraint(m, x[1]*x[2]*x[3]*y[1] >= 10)
	@NLconstraint(m, x[1]*y[1]*y[2] >= 10)
	@NLconstraint(m, x[1]*y[1]*y[2]*y[3] >= 10)
	@NLconstraint(m, x[1]*x[2]*x[3]*y[1]*y[2]*y[3] >= 10)
	@NLconstraint(m, x[1]*x[2]*x[3]*x[4]*x[5] >= 1)
	@NLconstraint(m, y[1]*y[2]*y[4]*y[5]*y[1] >= 99)

	@objective(m, Min, x[1]+x[2]+x[3])

	return m
end

function bmpl_linearlifting(;solver=nothing)

	m = Model(solver=solver)

	@variable(m, x[1:5], Bin)
	@variable(m, y[1:5]>=0)

	@NLconstraint(m, 15*x[1]*(2*y[1]+y[2]) >= 10)
	@NLconstraint(m, 15*(x[1]+x[2]*y[1])*(y[1]+y[2]*x[1]*x[2]) + 10 >= 10)
	@NLconstraint(m, 15*x[2]*(y[1]+y[2]*x[3]*y[5]*x[4]) + 50 >= 10)

	@NLobjective(m, Min, x[1] * (y[1]*y[2]*y[3] + x[2]*x[3]))

	return m
end

function bpml_lnl(solver=nothing)

	m = Model(solver=solver)

	srand(10)
	@variable(m, X[1:5], Bin)
	@variable(m, 0.1<=Y[1:5]<=0.1+10*rand())
	@constraint(m,  sum(X) >= 3)
	@NLobjective(m, Min, sum(X[i]*Y[i] for i in 1:5))

	return m
end

function bpml_binl(solver=nothing)

	m = Model(solver=solver)

	srand(10)
	@variable(m, X[1:5], Bin)
	@variable(m, 50<=Y[1:5]<=50+100*rand()*rand())
	@constraint(m, sum(X) >= 3)
	@NLobjective(m, Max, sum(X[i]*Y[i]*Y[i+1] for i in 1:4) - X[5]*Y[5]*Y[1])

	return m
end

function bpml_monl(solver=nothing)

	m = Model(solver=solver)

	srand(10)
	@variable(m, X[1:5], Bin)
	@variable(m, 50<=Y[1:5]<=50+100*rand()*rand())
	@constraint(m, sum(X) >= 3)
	@NLobjective(m, Max, sum(X[i]*Y[i]*Y[i] for i in 1:5))

	return m
end
