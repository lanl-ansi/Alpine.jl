function circle(;solver=nothing)

	m = Model(solver)

    @variable(m, 0<=x[1:2]<=2)
    @NLconstraint(m, x[1]^2 + x[2]^2 >= 2)
    @objective(m, Min, x[1]+x[2])

	return m
end

function circleN(;solver=nothing, N=2)

	m = Model(solver)

    @variable(m, 0<=x[1:N]<=N)
    @NLconstraint(m, sum(x[i]^2 for i in 1:N) >= N)
    @objective(m, Min, sum(x))

	return m
end
