@testset "Expression Operator Tests" begin
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

        @test m.internalModel.lifted_constr_expr_mip[1] == :(x[5]-1.0>=0.0)
        @test m.internalModel.lifted_constr_expr_mip[2] == :(x[6]-1.0<=0.0)
        @test m.internalModel.lifted_constr_expr_mip[3] == :(x[8]-1.0<=0.0)
        @test m.internalModel.lifted_constr_expr_mip[4] == :(x[10]-1.0>=0.0)
        @test m.internalModel.lifted_constr_expr_mip[5] == :(x[13]-1.0<=0.0)

        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[1])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[3])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[5]), :(x[7])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[9])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[3]), :(x[3])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_info, [:(x[5]), :(x[12])])

        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[1])]][:id] == 5
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[2])]][:id] == 6
        @test m.internalModel.nonlinear_info[[:(x[2]), :(x[2])]][:id] == 7
        @test m.internalModel.nonlinear_info[[:(x[5]), :(x[7])]][:id] == 8
        @test m.internalModel.nonlinear_info[[:(x[2]), :(x[3])]][:id] == 9
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[9])]][:id] == 10
        @test m.internalModel.nonlinear_info[[:(x[3]), :(x[3])]][:id] == 11
        @test m.internalModel.nonlinear_info[[:(x[7]), :(x[11])]][:id] == 12
        @test m.internalModel.nonlinear_info[[:(x[5]), :(x[12])]][:id] == 13
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

        @test m.internalModel.lifted_constr_expr_mip[1] == :(x[6]-1.0>=0.0)
        @test m.internalModel.lifted_constr_expr_mip[2] == :(x[11]-1.0<=0.0)
        @test m.internalModel.lifted_constr_expr_mip[3] == :(x[13]-1.0>=0.0)
        @test m.internalModel.lifted_constr_expr_mip[4] == :(x[15]-1.0<=0.0)
        @test m.internalModel.lifted_constr_expr_mip[5] == :(x[17]-1.0>=0.0)
        @test m.internalModel.lifted_constr_expr_mip[6] == :(x[19]-1.0<=0.0)

        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[2])]) #5
        @test haskey(m.internalModel.nonlinear_info, [:(x[5]), :(x[3])]) #6
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[1])]) #7
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[2])]) #8
        @test haskey(m.internalModel.nonlinear_info, [:(x[7]), :(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_info, [:(x[3]), :(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_info, [:(x[9]), :(x[10])]) #11
        @test haskey(m.internalModel.nonlinear_info, [:(x[8]), :(x[10])]) #12
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[12])]) #13
        @test haskey(m.internalModel.nonlinear_info, [:(x[7]), :(x[2])])  #14
        @test haskey(m.internalModel.nonlinear_info, [:(x[14]), :(x[10])]) #15
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[3])]) #16
        @test haskey(m.internalModel.nonlinear_info, [:(x[7]), :(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[8])]) #18
        @test haskey(m.internalModel.nonlinear_info, [:(x[18]), :(x[3])]) #19

        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[2])]][:id] == 5
        @test m.internalModel.nonlinear_info[[:(x[5]), :(x[3])]][:id] == 6
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[1])]][:id] == 7
        @test m.internalModel.nonlinear_info[[:(x[2]), :(x[2])]][:id] == 8
        @test m.internalModel.nonlinear_info[[:(x[7]), :(x[8])]][:id] == 9
        @test m.internalModel.nonlinear_info[[:(x[3]), :(x[3])]][:id] == 10
        @test m.internalModel.nonlinear_info[[:(x[9]), :(x[10])]][:id] == 11
        @test m.internalModel.nonlinear_info[[:(x[8]), :(x[10])]][:id] == 12
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[12])]][:id] == 13
        @test m.internalModel.nonlinear_info[[:(x[7]), :(x[2])]][:id] == 14
        @test m.internalModel.nonlinear_info[[:(x[14]), :(x[10])]][:id] == 15
        @test m.internalModel.nonlinear_info[[:(x[2]), :(x[3])]][:id] == 16
        @test m.internalModel.nonlinear_info[[:(x[7]), :(x[16])]][:id] == 17
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[8])]][:id] == 18
        @test m.internalModel.nonlinear_info[[:(x[18]), :(x[3])]][:id] == 19

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

        @test m.internalModel.lifted_constr_expr_mip[1] == :(x[7]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[2] == :(x[11]-1.0 <= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[3] == :(x[15]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[4] == :(x[18]-1.0 <= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[5] == :(x[20]-1.0 >= 0.0)
        @test m.internalModel.lifted_constr_expr_mip[6] == :(x[23]-1.0 <= 0.0)

        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[2])]) #5
        @test haskey(m.internalModel.nonlinear_info, [:(x[5]), :(x[3])]) #6
        @test haskey(m.internalModel.nonlinear_info, [:(x[6]), :(x[4])]) #7
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[1])]) #8
        @test haskey(m.internalModel.nonlinear_info, [:(x[8]), :(x[2])]) #9
        @test haskey(m.internalModel.nonlinear_info, [:(x[9]), :(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_info, [:(x[10]), :(x[4])]) #11
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[2])]) #12
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[12])])  #13
        @test haskey(m.internalModel.nonlinear_info, [:(x[13]), :(x[3])]) #14
        @test haskey(m.internalModel.nonlinear_info, [:(x[14]), :(x[4])]) #15
        @test haskey(m.internalModel.nonlinear_info, [:(x[3]), :(x[3])]) #16
        @test haskey(m.internalModel.nonlinear_info, [:(x[5]), :(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_info, [:(x[17]), :(x[4])]) #18
        @test haskey(m.internalModel.nonlinear_info, [:(x[4]), :(x[4])]) #19
        @test haskey(m.internalModel.nonlinear_info, [:(x[6]), :(x[19])]) #20
        @test haskey(m.internalModel.nonlinear_info, [:(x[8]), :(x[12])]) #21
        @test haskey(m.internalModel.nonlinear_info, [:(x[21]), :(x[16])]) #22
        @test haskey(m.internalModel.nonlinear_info, [:(x[22]), :(x[19])]) #23

        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[2])]][:id] == 5
        @test m.internalModel.nonlinear_info[[:(x[5]), :(x[3])]][:id] == 6
        @test m.internalModel.nonlinear_info[[:(x[6]), :(x[4])]][:id] == 7
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[1])]][:id] == 8
        @test m.internalModel.nonlinear_info[[:(x[8]), :(x[2])]][:id] == 9
        @test m.internalModel.nonlinear_info[[:(x[9]), :(x[3])]][:id] == 10
        @test m.internalModel.nonlinear_info[[:(x[10]), :(x[4])]][:id] == 11
        @test m.internalModel.nonlinear_info[[:(x[2]), :(x[2])]][:id] == 12
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[12])]][:id] == 13
        @test m.internalModel.nonlinear_info[[:(x[13]), :(x[3])]][:id] == 14
        @test m.internalModel.nonlinear_info[[:(x[14]), :(x[4])]][:id] == 15
        @test m.internalModel.nonlinear_info[[:(x[3]), :(x[3])]][:id] == 16
        @test m.internalModel.nonlinear_info[[:(x[5]), :(x[16])]][:id] == 17
        @test m.internalModel.nonlinear_info[[:(x[17]), :(x[4])]][:id] == 18
        @test m.internalModel.nonlinear_info[[:(x[4]), :(x[4])]][:id] == 19
        @test m.internalModel.nonlinear_info[[:(x[6]), :(x[19])]][:id] == 20
        @test m.internalModel.nonlinear_info[[:(x[8]), :(x[12])]][:id] == 21
        @test m.internalModel.nonlinear_info[[:(x[21]), :(x[16])]][:id] == 22
        @test m.internalModel.nonlinear_info[[:(x[22]), :(x[19])]][:id] == 23

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

        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[2]),:(x[3]),:(x[4])])) #5
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[1])]) #6
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[6]),:(x[2]),:(x[3]),:(x[4])])) #7
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[2])]) #8
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[8]),:(x[3]),:(x[4])])) #9
        @test haskey(m.internalModel.nonlinear_info, [:(x[3]), :(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_info, [:(x[4]), :(x[4])]) #11
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[2]),:(x[10]),:(x[11])])) #12
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[8]),:(x[10]),:(x[4])])) #13
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[6]),:(x[2]),:(x[3]),:(x[11])])) #14
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[6]),:(x[8]),:(x[10]),:(x[11])])) #15

        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]),:(x[2]),:(x[3]),:(x[4])])][:id] == 5 #5
        @test m.internalModel.nonlinear_info[[:(x[1]), :(x[1])]][:id] == 6 #6
        @test m.internalModel.nonlinear_info[Set(Any[:(x[6]),:(x[2]),:(x[3]),:(x[4])])][:id] == 7 #7
        @test m.internalModel.nonlinear_info[[:(x[2]), :(x[2])]][:id] == 8 #8
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]),:(x[8]),:(x[3]),:(x[4])])][:id] == 9 #9
        @test m.internalModel.nonlinear_info[[:(x[3]), :(x[3])]][:id] == 10 #10
        @test m.internalModel.nonlinear_info[[:(x[4]), :(x[4])]][:id] == 11 #11
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]),:(x[2]),:(x[10]),:(x[11])])][:id] ==  12  #12
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]),:(x[8]),:(x[10]),:(x[4])])][:id] == 13 #13
        @test m.internalModel.nonlinear_info[Set(Any[:(x[6]),:(x[2]),:(x[3]),:(x[11])])][:id] == 14 #14
        @test m.internalModel.nonlinear_info[Set(Any[:(x[6]),:(x[8]),:(x[10]),:(x[11])])][:id] == 15 #15

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

        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[2]),:(x[3])])) #5
        @test haskey(m.internalModel.nonlinear_info, [:(x[5]), :(x[4])]) #6
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[1])]) #7
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[7]),:(x[2]),:(x[3])])) #8
        @test haskey(m.internalModel.nonlinear_info, [:(x[8]), :(x[4])]) #9
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[2])]) #10
        @test haskey(m.internalModel.nonlinear_info, [:(x[10]), :(x[3])]) #11
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[11]),:(x[4])])) #12
        @test haskey(m.internalModel.nonlinear_info, [:(x[3]), :(x[3])]) #13
        @test haskey(m.internalModel.nonlinear_info, [:(x[2]), :(x[13])]) #14
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[1]),:(x[14]),:(x[4])])) #15
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[10])]) #16
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[16]),:(x[3]),:(x[4])])) #17
        @test haskey(m.internalModel.nonlinear_info, [:(x[1]), :(x[2])]) #18
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[18]),:(x[13]),:(x[4])])) #19
        @test haskey(m.internalModel.nonlinear_info, [:(x[4]), :(x[4])]) #20
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[18]),:(x[3]),:(x[20])])) #21
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[18]),:(x[13]),:(x[20])])) #22
        @test haskey(m.internalModel.nonlinear_info, Set(Any[:(x[7]),:(x[10]),:(x[13])])) #23
        @test haskey(m.internalModel.nonlinear_info, [:(x[23]), :(x[20])]) #24

    end

end
