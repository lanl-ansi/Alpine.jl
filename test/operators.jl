o@testset "Expression Operator Tests" begin
    @testset " Validation Test on expression parsing : part1 " begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(OutputFlag=0),
								   log_level=0))
        @variable(m, x[1:4]>=0)
        @NLconstraint(m, x[1]^2 >= 1)  					# Basic monomial x[5]=x[1]^2
        @NLconstraint(m, x[1]*x[2] <= 1)				# x[6] <= 1 : x[6] = x[1]*x[2]
        @NLconstraint(m, x[1]^2 * x[2]^2 <= 1)          # x[5] + x[7] <= 1 : x[7] = x[2]^2
        @NLconstraint(m, x[1]*(x[2]*x[3]) >= 1)         # x[9] >= 1 : x[8] = x[2] * x[3] && x[9] = x[1]*x[8]
        @NLconstraint(m, x[1]^2*(x[2]^2 * x[3]^2) <= 1) # x[12] <= 1 : x[10] = x[3] ^ 2 && x[11] = x[7] * x[10] && x[12] = x[11]*[5]

        JuMP.build(m)

        @test m.internalModel.bounding_constr_expr_mip[1] == :(x[5]-1.0>=0.0)
        @test m.internalModel.bounding_constr_expr_mip[2] == :(x[6]-1.0<=0.0)
        @test m.internalModel.bounding_constr_expr_mip[3] == :(x[8]-1.0<=0.0)
        @test m.internalModel.bounding_constr_expr_mip[4] == :(x[10]-1.0>=0.0)
        @test m.internalModel.bounding_constr_expr_mip[5] == :(x[13]-1.0<=0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[1])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[3])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]), :(x[7])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[9])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[3])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]), :(x[12])])

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:id] == 3
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[7])]][:id] == 4
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[3])]][:id] == 5
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[9])]][:id] == 6
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[3])]][:id] == 7
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[11])]][:id] == 8
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[12])]][:id] == 9
    end

    @testset " Validation Test on expression parsing : part2" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0))

        @variable(m, x[1:4]>=0)
        @NLconstraint(m, (x[1]*x[2]) * x[3] >= 1)
        @NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

        @NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
        @NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2 <= 1)
        @NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
        @NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

        JuMP.build(m)

        @test m.internalModel.bounding_constr_expr_mip[1] == :(x[6]-1.0>=0.0)
        @test m.internalModel.bounding_constr_expr_mip[2] == :(x[11]-1.0<=0.0)
        @test m.internalModel.bounding_constr_expr_mip[3] == :(x[13]-1.0>=0.0)
        @test m.internalModel.bounding_constr_expr_mip[4] == :(x[15]-1.0<=0.0)
        @test m.internalModel.bounding_constr_expr_mip[5] == :(x[17]-1.0>=0.0)
        @test m.internalModel.bounding_constr_expr_mip[6] == :(x[19]-1.0<=0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[2])]) #5
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]), :(x[3])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[1])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[2])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]), :(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[9]), :(x[10])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[8]), :(x[10])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[12])]) #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]), :(x[2])])  #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[14]), :(x[10])]) #15
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[3])]) #16
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]), :(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[8])]) #18
        @test haskey(m.internalModel.nonlinear_terms, [:(x[18]), :(x[3])]) #19

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[3])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:id] == 3
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:id] == 4
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[8])]][:id] == 5
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[3])]][:id] == 6
        @test m.internalModel.nonlinear_terms[[:(x[9]), :(x[10])]][:id] == 7
        @test m.internalModel.nonlinear_terms[[:(x[8]), :(x[10])]][:id] == 8
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[12])]][:id] == 9
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[2])]][:id] == 10
        @test m.internalModel.nonlinear_terms[[:(x[14]), :(x[10])]][:id] == 11
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[3])]][:id] == 12
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[16])]][:id] == 13
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[8])]][:id] == 14
        @test m.internalModel.nonlinear_terms[[:(x[18]), :(x[3])]][:id] == 15
    end

    @testset "Validation Test on expression parsing: part3" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0))

        @variable(m, x[1:4]>=0)
        @NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 1)
        @NLconstraint(m, ((x[1]^2*x[2])*x[3])*x[4] <= 1)
        @NLconstraint(m, ((x[1]*x[2]^2)*x[3])*x[4] >= 1)
        @NLconstraint(m, ((x[1]*x[2])*x[3]^2)*x[4] <= 1)
        @NLconstraint(m, ((x[1]*x[2])*x[3])*x[4]^2 >= 1)
        @NLconstraint(m, ((x[1]^2*x[2]^2)*x[3]^2)*x[4]^2 <= 1)

        JuMP.build(m)

        @test m.internalModel.bounding_constr_expr_mip[1] == :(x[7]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[2] == :(x[11]-1.0 <= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[3] == :(x[15]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[4] == :(x[18]-1.0 <= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[5] == :(x[20]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[6] == :(x[23]-1.0 <= 0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[2])]) #5
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]), :(x[3])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]), :(x[4])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[1])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[8]), :(x[2])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[9]), :(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[10]), :(x[4])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[2])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[12])])  #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[13]), :(x[3])]) #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[14]), :(x[4])]) #15
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[3])]) #16
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]), :(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_terms, [:(x[17]), :(x[4])]) #18
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[4])]) #19
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]), :(x[19])]) #20
        @test haskey(m.internalModel.nonlinear_terms, [:(x[8]), :(x[12])]) #21
        @test haskey(m.internalModel.nonlinear_terms, [:(x[21]), :(x[16])]) #22
        @test haskey(m.internalModel.nonlinear_terms, [:(x[22]), :(x[19])]) #23

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[3])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[6]), :(x[4])]][:id] == 3
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:id] == 4
        @test m.internalModel.nonlinear_terms[[:(x[8]), :(x[2])]][:id] == 5
        @test m.internalModel.nonlinear_terms[[:(x[9]), :(x[3])]][:id] == 6
        @test m.internalModel.nonlinear_terms[[:(x[10]), :(x[4])]][:id] == 7
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:id] == 8
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[12])]][:id] == 9
        @test m.internalModel.nonlinear_terms[[:(x[13]), :(x[3])]][:id] == 10
        @test m.internalModel.nonlinear_terms[[:(x[14]), :(x[4])]][:id] == 11
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[3])]][:id] == 12
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[16])]][:id] == 13
        @test m.internalModel.nonlinear_terms[[:(x[17]), :(x[4])]][:id] == 14
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[4])]][:id] == 15
        @test m.internalModel.nonlinear_terms[[:(x[6]), :(x[19])]][:id] == 16
        @test m.internalModel.nonlinear_terms[[:(x[8]), :(x[12])]][:id] == 17
        @test m.internalModel.nonlinear_terms[[:(x[21]), :(x[16])]][:id] == 18
        @test m.internalModel.nonlinear_terms[[:(x[22]), :(x[19])]][:id] == 19
    end

    @testset "Validation Test on expression parsing: part7" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0))
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]*x[2]*x[3]^2*x[4]^2 >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]^2*x[4] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]*x[4]^2 >= 1)
        @NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4]^2 >= 1)

        JuMP.build(m)

        @test m.internalModel.lifted_constr_expr_mip[1] == :(x[5]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[2] == :(x[7]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[3] == :(x[9]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[4] == :(x[12]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[5] == :(x[13]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[6] == :(x[14]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[7] == :(x[15]-1.0 >= 0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2]),:(x[3]),:(x[4])]) #5
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[1])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]),:(x[2]),:(x[3]),:(x[4])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[2])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[3]),:(x[4]),:(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[4])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2]),:(x[10]),:(x[11])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[4]),:(x[8]),:(x[10])]) #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]),:(x[11]),:(x[2]),:(x[3])]) #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]),:(x[8]),:(x[10]),:(x[11])]) #15

        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[2]),:(x[3]),:(x[4])]][:id] == 1 #5
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[1])]][:id] == 2 #6
        @test m.internalModel.nonlinear_terms[[:(x[6]),:(x[2]),:(x[3]),:(x[4])]][:id] == 3 #7
        @test m.internalModel.nonlinear_terms[[:(x[2]),:(x[2])]][:id] == 4 #8
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[3]),:(x[4]),:(x[8])]][:id] == 5 #9
        @test m.internalModel.nonlinear_terms[[:(x[3]),:(x[3])]][:id] == 6 #10
        @test m.internalModel.nonlinear_terms[[:(x[4]),:(x[4])]][:id] == 7 #11
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[2]),:(x[10]),:(x[11])]][:id] ==  8  #12
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[4]),:(x[8]),:(x[10])]][:id] == 9 #13
        @test m.internalModel.nonlinear_terms[[:(x[6]),:(x[11]),:(x[2]),:(x[3])]][:id] == 10 #14
        @test m.internalModel.nonlinear_terms[[:(x[6]),:(x[8]),:(x[10]),:(x[11])]][:id] == 11 #15
    end

    @testset "Validation Test on expression parsing: part8" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0))
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, (x[1]*x[2]*x[3])*x[4] >= 1)
        @NLconstraint(m, (x[1]^2*x[2]*x[3])*x[4] >= 1)
        @NLconstraint(m, x[1]*(x[2]^2*x[3])*x[4] >= 1)
        @NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2]^2)*x[3]*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2])*x[3]*x[4]^2 >= 1)
        @NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4]^2 >= 1)
        @NLconstraint(m, (x[1]^2*x[2]^2*x[3]^2)*x[4]^2 >= 1)

        JuMP.build(m)

        @test m.internalModel.lifted_constr_expr_mip[1] == :(x[6]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[2] == :(x[9]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[3] == :(x[12]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[4] == :(x[15]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[5] == :(x[17]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[6] == :(x[19]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[7] == :(x[21]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[8] == :(x[22]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[9] == :(x[24]-1.0 >= 0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2]),:(x[3])]) #5
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]),:(x[4])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[1])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]),:(x[2]),:(x[3])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[8]),:(x[4])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[2])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[10]),:(x[3])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[4]),:(x[11])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[3])]) #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[13])]) #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[4]),:(x[14])]) #15
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[10])]) #16
        @test haskey(m.internalModel.nonlinear_terms, [:(x[16]),:(x[3]),:(x[4])]) #17
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2])]) #18
        @test haskey(m.internalModel.nonlinear_terms, [:(x[18]),:(x[13]),:(x[4])]) #19
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[4])]) #20
        @test haskey(m.internalModel.nonlinear_terms, [:(x[18]),:(x[20]),:(x[3])]) #21
        @test haskey(m.internalModel.nonlinear_terms, [:(x[18]),:(x[13]),:(x[20])]) #22
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]),:(x[10]),:(x[13])]) #23
        @test haskey(m.internalModel.nonlinear_terms, [:(x[23]),:(x[20])]) #24
    end

    @testset "Validation Test on constraint parsing: part1" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                           mip_solver=CbcSolver(OutputFlag=0),
                           log_level=0))


        @variable(m, 0<=x[1:5]<=2)

        @constraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                             # 1: true
        @constraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                         # 2: true
        @constraint(m, 3(x[1]x[1]) + 4*x[2]*x[2] <= -5)                             # 3: false
        @constraint(m, 3(x[1]x[1]) + 4*x[2]^2 <= 10)                                # 4: true
        @constraint(m, 3x[1]^2 + 4x[2]^2 + 6x[3]^2 <= 10)                           # 5: true

        @NLconstraint(m, 3x[1]^0.5 + 4x[2]^0.5 + 5x[5]^0.5 <= 100)					# 6: true | type-C
        @NLconstraint(m, -3x[1]^0.5 -4x[2]^0.5 >= -100)								# 7: true | type-C
        @NLconstraint(m, 3x[1]^3 + x[2]^3 + 5x[3]^3 <= 200)                         # 8: true | type-a
        @NLconstraint(m, x[1]*x[1]*x[1] + x[2]*x[2]*x[2] + x[3]*x[3]*x[3] + 5*x[4]*20*x[4]*x[4] <= 200) # 9: true
        @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                           # 10: true

        @NLconstraint(m, (3*x[1]*x[1] + 4*x[2]*x[2]) <= 25)                         # 11: true
        @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] - 25 <= 0)                       # 12: true
        @NLconstraint(m, -3*x[1]*x[1] -4*x[2]*x[2] >= -25)                          # 13: true
        @NLconstraint(m, 3*x[1]*x[1] + 5x[2]*x[2] <= 25)                            # 14: true
        @NLconstraint(m, x[1]*3*x[1] + x[2]*x[2]*5 + x[4]^(3-1) <= 25)              # 15: true

        @NLconstraint(m, 4*x[1]^2 + 5x[2]^2 <= 25)                                  # 16: true
        @NLconstraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                       # 17: false (unsupported when with @NLconstraint)
        @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[1] <= 25)                           # 18: false
        @NLconstraint(m, 3*x[1]*x[1] + 16*x[2]^2 <= 40)                             # 19: true
        @NLconstraint(m, 3*x[1]^2 + 16*x[2]^2 + 17 <= 16)                           # 20: false

        @NLconstraint(m, 3*x[1]^3 + 16*x[2]^2 <= 20 - 20)                           # 21: false
        @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] + 5*x[3]*x[3] + 6x[4]x[4] <= 15) # 22: true
        @NLconstraint(m, 3x[1]x[1] + 4x[2]x[2] + 5x[3]^2 <= -15)                    # 23: false
        @NLconstraint(m, 3x[1]^2 + 4x[2]^2 >= 15)                                   # 24: false
        @NLconstraint(m, - 3x[1]^2 - 4x[3]^2 >= -15)                                # 25: true

        @NLconstraint(m, 3x[1]^4 + 4x[2]^4 <= 200)                                  # 26: true
        @NLconstraint(m, 3x[1]^4 + 4x[2]x[2]x[2]x[2] - 200 <= 0)                    # 27: true
        @NLconstraint(m, 3x[1]^4 + 4x[2]^2*x[2]*x[2] <= 200)                        # 28: true
        @NLconstraint(m, 3x[1]^4 + 4x[2]^3 <= 200)                                  # 29: false
        @NLconstraint(m, 3x[1]^8 + 16*25*x[2]^8 - 30x[3]^8 <= 50)                   # 30: false

        @objective(m, Max, x[1]^2+x[3]^2)                                           # true

        JuMP.build(m)

        @test m.internalModel.num_constr_convex == 21
        @test length(m.internalModel.nonlinear_terms) ==

        # 0 : OBJ
        @test m.internalModel.structural_obj == :convex
        @test m.internalModel.nonlinear_constrs[0][:expr_orig] == :objecive
        @test m.internalModel.nonlinear_constrs[0][:idx] == 0
        @test m.internalModel.nonlinear_constrs[0][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[0][:convexified] == :false

        @test m.internalModel.bounding_obj_mip[:sense] == nothing
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0, 1.0]
        @test m.internalModel.bounding_obj_mip[:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_obj_mip[:rhs] == 0.0
        @test m.internalModel.bounding_obj_mip[:powers] == [2, 2]
        @test m.internalModel.bounding_obj_mip[:cnt] == 2

        # 1
        @test m.internalModel.structural_constr[1] == :convex
        @test m.internalModel.nonlinear_constrs[1][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[1][:idx] == 0
        @test m.internalModel.nonlinear_constrs[1][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[1][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[:coefs] == [3.0, 4.0]
        @test m.internalModel.bounding_constr_mip[:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[:cnt] == 2

        # 2...
    end

end
