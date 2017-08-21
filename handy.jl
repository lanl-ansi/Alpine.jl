using JuMP, MathProgBase, POD, Cbc, Ipopt

m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                   mip_solver=CbcSolver(OutputFlag=0),
                   log_level=0))

@variable(m, 0<=x[1:5]<=2)

@constraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                             # 1: :convex
@constraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                         # 2: :convex
@constraint(m, 3(x[1]x[1]) + 4*x[2]*x[2] <= -5)                             # 3: :affine
@constraint(m, 3(x[1]x[1]) + 4*x[2]^2 <= 10)                                # 4: :convex
@constraint(m, 3x[1]^2 + 4x[2]^2 + 6x[3]^2 <= 10)                           # 5: :convex

@NLconstraint(m, 3x[1]^0.5 + 4x[2]^0.5 + 5x[5]^0.5 <= 100)					# 6: :convex | type-C
@NLconstraint(m, -3x[1]^0.5 -4x[2]^0.5 >= -100)								# 7: :convex | type-C
@NLconstraint(m, 3x[1]^3 + x[2]^3 + 5x[3]^3 <= 200)                         # 8: :convex | type-a
@NLconstraint(m, x[1]*x[1]*x[1] + x[2]*x[2]*x[2] + x[3]*x[3]*x[3] + 5*x[4]*20*x[4]*x[4] <= 200) # 9: :convex
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                           # 10: :convex

@NLconstraint(m, (3*x[1]*x[1] + 4*x[2]*x[2]) <= 25)                         # 11: :convex
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] - 25 <= 0)                       # 12: :convex
@NLconstraint(m, -3*x[1]*x[1] -4*x[2]*x[2] >= -25)                          # 13: :convex
@NLconstraint(m, 3*x[1]*x[1] + 5x[2]*x[2] <= 25)                            # 14: :convex
@NLconstraint(m, x[1]*3*x[1] + x[2]*x[2]*5 + x[4]^(3-1) <= 25)              # 15: :convex

@NLconstraint(m, 4*x[1]^2 + 5x[2]^2 <= 25)                                  # 16: :convex
@NLconstraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                       # 17: :affine (unsupported when with @NLconstraint)
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[1] <= 25)                           # 18: :affine
@NLconstraint(m, 3*x[1]*x[1] + 16*x[2]^2 <= 40)                             # 19: :convex
@NLconstraint(m, 3*x[1]^2 + 16*x[2]^2 + 17 <= 16)                           # 20: :affine

@NLconstraint(m, 3*x[1]^3 + 16*x[2]^2 <= 20 - 20)                           # 21: :affine
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] + 5*x[3]*x[3] + 6x[4]x[4] <= 15) # 22: :convex
@NLconstraint(m, 3x[1]x[1] + 4x[2]x[2] + 5x[3]^2 <= -15)                    # 23: :affine
@NLconstraint(m, 3x[1]^2 + 4x[2]^2 >= 15)                                   # 24: :affine
@NLconstraint(m, - 3x[1]^2 - 4x[3]^2 >= -15)                                # 25: :convex

@NLconstraint(m, 3x[1]^4 + 4x[2]^4 <= 200)                                  # 26: :convex
@NLconstraint(m, 3x[1]^4 + 4x[2]x[2]x[2]x[2] - 200 <= 0)                    # 27: :convex
@NLconstraint(m, 3x[1]^4 + 4x[2]^2*x[2]*x[2] <= 200)                        # 28: :convex
@NLconstraint(m, 3x[1]^4 + 4x[2]^3 <= 200)                                  # 29: :affine
@NLconstraint(m, 3x[1]^8 + 16*25*x[2]^8 - 30x[3]^8 <= 50)                   # 30: :affine

@objective(m, Max, x[1]^2+x[3]^2)

JuMP.build(m)

test_range = 30

test_map = Dict(true=>"pass", false=>"fail")

test_answer = [:convex, :convex, :affine, :convex, :convex,
               :convex, :convex, :convex, :convex, :convex,
               :convex, :convex, :convex, :convex, :convex,
               :convex, :affine, :affine, :convex, :affine,
               :affine, :convex, :affine, :affine, :convex,
               :convex, :convex, :convex, :affine, :affine]

println("Test | Convex? | Expression")
for ex_i in 1:test_range
    println("$(test_map[m.internalModel.structural_constr[ex_i] == test_answer[ex_i]]) | $(m.internalModel.bounding_constr_expr_mip[ex_i])")
end
