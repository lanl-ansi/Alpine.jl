using JuMP, MathProgBase, POD

m = Model()

@variable(m, x[1:5])
@constraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                             # true
@constraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                         # true
@constraint(m, 3(x[1]x[1]) + 4*x[2]*x[2] <= -5)                             # false
@constraint(m, 3(x[1]x[1]) + 4*x[2]^2 <= 10)                                # true
@constraint(m, 3x[1]^2 + 4x[2]^2 + 6x[3]^2 <= 10)                           # true

@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                           # true
@NLconstraint(m, (3*x[1]*x[1] + 4*x[2]*x[2]) <= 25)                         # true
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] - 25 <= 0)                       # true
@NLconstraint(m, -3*x[1]*x[1] -4*x[2]*x[2] >= -25)                          # true
@NLconstraint(m, 3*x[1]*x[1] + 5x[2]*x[2] <= 25)                            # true

@NLconstraint(m, 4*x[1]^2 + 5x[2]^2 <= 25)                                  # Pass
@NLconstraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                       # false (unsupported when with @NLconstraint)
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[1] <= 25)                           # false
@NLconstraint(m, 3*x[1]*x[1] + 16*x[2]^2 <= 40)                             # true
@NLconstraint(m, 3*x[1]^2 + 16*x[2]^2 + 17 <= 16)                           # false

@NLconstraint(m, 3*x[1]^3 + 16*x[2]^2 <= 20 - 20)                           # false
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] + 5*x[3]*x[3] + 6x[4]x[4] <= 15) # true
@NLconstraint(m, 3x[1]x[1] + 4x[2]x[2] + 5x[3]^2 <= -15)                    # false
@NLconstraint(m, 3x[1]^2 + 4x[2]^2 >= 15)                                   # false
@NLconstraint(m, - 3x[1]^2 - 4x[3]^2 >= -15)                                # true

@NLconstraint(m, 3x[1]^4 + 4x[2]^4 <= 200)                                  # true
@NLconstraint(m, 3x[1]^4 + 4x[2]x[2]x[2]x[2] - 200 <= 0)                    # true
@NLconstraint(m, 3x[1]^4 + 4x[2]^2*x[2]*x[2] <= 200)                        # true
@NLconstraint(m, 3x[1]^4 + 4x[2]^3 <= 200)                                  # false
@NLconstraint(m, 3x[1]^8 + 16*25*x[2]^8 - 30x[3]^8 <= 50)                   # false

@NLconstraint(m, )

test_range = 25

d = JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:ExprGraph])

exs = []
exs_convex = []
for ex_i in 1:test_range
    push!(exs, MathProgBase.constr_expr(d, ex_i))
    POD.preprocess_expression(exs[end])
    push!(exs_convex, POD.resolve_convex_constr(exs[end]))
    println("-------")
end

test_map = Dict(true=>"pass", false=>"fail")
test_answer = [true, true, false, true, true,
               true, true, true, true, true,
               true, false, false, true, false,
               false, true, false, false, true,
               true, true, true, false, false]

println("Test | Convex? | Expression")
for ex_i in 1:test_range
    if exs_convex[ex_i]
        println("$(test_map[exs_convex[ex_i] == test_answer[ex_i]]) | $(exs_convex[ex_i])    | $(exs[ex_i])")
    else
        println("$(test_map[exs_convex[ex_i] == test_answer[ex_i]]) | $(exs_convex[ex_i])   | $(exs[ex_i])")
    end
end
