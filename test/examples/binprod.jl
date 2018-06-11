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
