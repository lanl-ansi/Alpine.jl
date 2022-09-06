function basic_linear_lift(; verbose = false, solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, x[i = 1:3] >= 1) # At some point if an initial value is given, keep them
    @constraint(m, (x[1] - x[2]) * (3 * x[2] - x[3]) >= 111)
    @NLconstraint(m, (x[1] - x[2]) * (3 * x[2] - x[3]) >= 111)
    @NLconstraint(m, (x[1] - x[2])^2 >= 100)
    @NLconstraint(m, (x[1] - x[2]) * (3 * x[2] - x[3]) * (x[1] + x[3]) >= 111)
    @NLconstraint(m, (3 + x[1] + x[2] + x[3]) * (x[1] + x[2])^2 >= 111)

    # Next level goal
    # @NLconstraint(m, (9-x[1]-x[2]-x[3])^2 >= 111)

    @objective(m, Min, x[1] + x[2] + x[3])

    return m
end
