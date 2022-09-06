function div(; verbose = false, solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, 1 <= x[1:2] <= 10)

    @NLobjective(m, Min, 6 * x[1]^2 + 4 * x[2]^2 - 2.5 / 5 * x[1] * x[2])

    @constraint(m, x[1] * (-1 * 3 / 5 * 5) >= 0)
    @constraint(m, x[2] * 5 / 20 >= 0)
    @constraint(m, x[2] * (5 * 15 / 15) >= 0)
    @constraint(m, x[2] * (120 * 2 / -2) >= 0)
    @constraint(m, 1 * 2 * 3 * 4 * 5 * 6 * x[2] / 0.01 >= 0)
    @constraint(m, x[1] * 1 * 2 * 3 * 4 * 5 * 6 / 0.01 >= 0)
    @NLconstraint(m, (3 / 5) * x[1] * (60 / 60) * x[2] >= 8)
    @NLconstraint(
        m,
        (1 * 2 * 3 * 4 / 5 / 6 * 7) * x[2] - x[1] * 1 * 2 * 3 * 4 * 5 * 6 * x[2] / 0.01 >=
        0
    )
    @NLconstraint(
        m,
        (1 * 2 * 3 * 4 / 5 / 6 * 7) * x[2] -
        0.5(x[1] * 1 * 2 * 3 * 4 * 5 * 6 * x[2] / 0.01) >= 0
    )
    @NLconstraint(
        m,
        (1 * 2 * 3 * 4 / 5 / 6 * 7) * x[2] -
        0.5(x[1] * 2 * 3 * x[2] / 0.01 + 5 * 7 / 10 * x[2]) >= 0
    )
    @NLconstraint(m, (1 * 2 * 3 * 4 / 5 / 6 * 7) * x[2] - 0.5 * (x[1] - x[2]x[1]) >= 0)

    return m
end
