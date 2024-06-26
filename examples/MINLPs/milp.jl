function milp(; solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, 0 <= objvar <= 20, Int)
    @objective(m, Min, objvar)

    return m
end
