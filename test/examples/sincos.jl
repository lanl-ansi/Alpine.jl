function sincos_p1(;solver=nothing)

    m = Model(solver)

    @variable(m, 0.1<=x[1:10]<=10)
    @NLconstraint(m, [i in 1:9], x[i]*x[i+1]*sin(x[i])>=0.32)
    @NLconstraint(m, [i in 1:9], x[i]*x[i+1]*sin(x[i])<=0.42)
    @NLobjective(m, Min, sum(cos(x[i])+sin(x[i]) for i in 1:10))

    return m
end

function trig(;solver=nothing)

    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1]
    @variable(m, x[x_Idx])
    JuMP.set_lower_bound(x[1], -2.0)
    JuMP.set_upper_bound(x[1], 5.0)

    # ----- Constraints ----- #
    @NLconstraint(m, e1, -(sin(11*x[1])+cos(13*x[1])-sin(17*x[1])-cos(19*x[1]))+objvar == 0.0)
    @NLconstraint(m, e2, 5*sin(x[1])-x[1] <= 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function specialopts(;verbose=false, solver=nothing)

	m = Model(solver)

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

	return m
end
