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
