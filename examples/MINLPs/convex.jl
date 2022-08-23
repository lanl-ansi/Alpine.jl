function convex_test(; solver = nothing)
    m = Model(solver)

    @variable(m, 0 <= x[1:5] <= 2)

    @constraint(m, 3 * x[1] * x[1] + 4 * x[2] * x[2] <= 25)                             # 1: true
    @constraint(m, 3 * x[1] * x[1] - 25 + 4 * x[2] * x[2] <= 0)                         # 2: true
    @constraint(m, 3(x[1]x[1]) + 4 * x[2] * x[2] <= -5)                             # 3: false
    @constraint(m, 3(x[1]x[1]) + 4 * x[2]^2 <= 10)                                # 4: true
    @constraint(m, 3x[1]^2 + 4x[2]^2 + 6x[3]^2 <= 10)                           # 5: true

    @NLconstraint(m, 3x[1]^0.5 + 4x[2]^0.5 + 5x[5]^0.5 <= 100)# 6: true | type-C
    @NLconstraint(m, -3x[1]^0.5 - 4x[2]^0.5 >= -100)# 7: true | type-C
    @NLconstraint(m, 3x[1]^3 + x[2]^3 + 5x[3]^3 <= 200)                         # 8: true | type-a
    @NLconstraint(
        m,
        x[1] * x[1] * x[1] +
        x[2] * x[2] * x[2] +
        x[3] * x[3] * x[3] +
        5 * x[4] * 20 * x[4] * x[4] <= 200
    ) # 9: true
    @NLconstraint(m, 3 * x[1] * x[1] + 4 * x[2] * x[2] <= 25)                           # 10: true

    @NLconstraint(m, (3 * x[1] * x[1] + 4 * x[2] * x[2]) <= 25)                         # 11: true
    @NLconstraint(m, 3 * x[1] * x[1] + 4 * x[2] * x[2] - 25 <= 0)                       # 12: true
    @NLconstraint(m, -3 * x[1] * x[1] - 4 * x[2] * x[2] >= -25)                          # 13: true
    @NLconstraint(m, 3 * x[1] * x[1] + 5x[2] * x[2] <= 25)                            # 14: true
    @NLconstraint(m, x[1] * 3 * x[1] + x[2] * x[2] * 5 + x[4]^(3 - 1) <= 25)              # 15: true

    @NLconstraint(m, 4 * x[1]^2 + 5x[2]^2 <= 25)                                  # 16: true
    @NLconstraint(m, 3 * x[1] * x[1] - 25 + 4 * x[2] * x[2] <= 0)                       # 17: false (unsupported when with @NLconstraint)
    @NLconstraint(m, 3 * x[1] * x[1] + 4 * x[2] * x[1] <= 25)                           # 18: false
    @NLconstraint(m, 3 * x[1] * x[1] + 16 * x[2]^2 <= 40)                             # 19: true
    @NLconstraint(m, 3 * x[1]^2 + 16 * x[2]^2 + 17 <= 16)                           # 20: false

    @NLconstraint(m, 3 * x[1]^3 + 16 * x[2]^2 <= 20 - 20)                           # 21: false
    @NLconstraint(
        m,
        3 * x[1] * x[1] + 4 * x[2] * x[2] + 5 * x[3] * x[3] + 6x[4]x[4] <= 15
    ) # 22: true
    @NLconstraint(m, 3x[1]x[1] + 4x[2]x[2] + 5x[3]^2 <= -15)                    # 23: false
    @NLconstraint(m, 3x[1]^2 + 4x[2]^2 >= 15)                                   # 24: false
    @NLconstraint(m, sum(x[i]^2 for i in 1:5) <= 99999)                         # 25: true

    @NLconstraint(m, 3x[1]^4 + 4x[2]^4 <= 200)                                  # 26: true
    @NLconstraint(m, 3x[1]^4 + 4x[2]x[2]x[2]x[2] - 200 <= 0)                    # 27: true
    @NLconstraint(m, 3x[1]^4 + 4x[2]^2 * x[2] * x[2] <= 200)                        # 28: true
    @NLconstraint(m, 3x[1]^4 + 4x[2]^3 <= 200)                                  # 29: false
    @NLconstraint(m, 3x[1]^8 + 16 * 25 * x[2]^8 - 30x[3]^8 <= 50)                   # 30: false

    @objective(m, Max, x[1]^2 + x[3]^2)                                           # true

    return m
end

function convex_solve(; solver = nothing)
    m = Model(solver)

    @variable(m, 0 <= x[1:5] <= 100)

    @NLconstraint(m, 3 * x[1] * x[1] + 4 * x[2] * x[2] <= 25)                             # 1: true
    @NLconstraint(m, 3(x[1]x[1]) + 4 * x[2]^2 <= 10)                                # 4: true
    @NLconstraint(m, 3x[1]^2 + 4x[2]^2 + 6x[3]^2 <= 10)                           # 5: true
    @NLconstraint(m, (x[1] + x[2])^2 <= 100)
    @NLconstraint(m, -3 * x[1] * x[1] - 4 * x[2] * x[2] >= -25)                            # 13: true

    @objective(m, Max, x[1]^2 + x[3]^2)                                             # true

    return m
end
