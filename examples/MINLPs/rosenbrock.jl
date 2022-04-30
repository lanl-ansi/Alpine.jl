# Reference: https://direct.mit.edu/evco/article/14/1/119/1232/A-Note-on-the-Extended-Rosenbrock-Function

function rosenbrock(;solver = nothing)

	m = Model(solver)

    @variable(m, -10 <= x <= 10)
    @variable(m, -10 <= y <= 10)

	@NLobjective(m, Min, 100*(y-x^2)^2 + (1-x)^2)

	return m
end