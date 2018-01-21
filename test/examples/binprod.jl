function circle(;verbose=false, solver=nothing)

	m = Model(solver=solver)

    @variable(m, x[1:5], Bin)
    @NLconstraint(m, x[1]^2 + x[2]^2 >= 2)
    @objective(m, Min, x[1]+x[2])

	if verbose
		print(m)
	end

	return m
end

function binprod_nlp3(;solver=nothing)

	m = Model(solver=solver)

	@variable(m, x[1:8])
	@variable(m, y[1:5], Bin)

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

	@constraint(m, 0.0025*(x[4]*y[1] + x[6]*y[2]) <= 1)
	@constraint(m, 0.0025*(x[5] - x[4]*y[1] + x[7]) <= 1)
	@constraint(m, 0.01(x[8]-x[5]*y[3]) <= 1)
	@NLconstraint(m, 100*x[1] - x[1]*x[6]*y[1] + 833.33252*x[4]*y[1] <= 83333.333)
	@NLconstraint(m, x[2]*x[4]*y[4] - x[2]*x[7] - 1250*x[4] + 1250*x[5] <= 0)
	@NLconstraint(m, x[3]*x[5]*y[2]*y[5] - x[3]*x[8]*y[5] - 2500*x[5]*y[1]*y[4] + 1250000 <= 0)

	@NLconstraint(m, y[1]*y[2]*y[3] <= 0)
	@NLconstraint(m, y[2]*y[3] >= y[4]*y[5])
	@NLconstraint(m, y[1]*y[5] <= y[2]*y[4])

	@objective(m, Min, x[1]+x[2]+x[3])

	return m
end
