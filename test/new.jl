@testset "Compact nlexpr.jl tests" begin

    # for i in 1:m.internalModel.num_constr_orig
    #     if m.internalModel.structural_constr[i] == :affine
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:rhs] == $(m.internalModel.bounding_constr_mip[i][:rhs])")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:vars] == $(m.internalModel.bounding_constr_mip[i][:vars])")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:coefs] == $(m.internalModel.bounding_constr_mip[i][:coefs])")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:sense] == :($(m.internalModel.bounding_constr_mip[i][:sense]))")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:cnt] == $(m.internalModel.bounding_constr_mip[i][:cnt])")
    #     elseif m.internalModel.structural_constr[i] == :convex
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:rhs] == $(m.internalModel.bounding_constr_mip[i][:rhs])")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:vars] == $(m.internalModel.bounding_constr_mip[i][:vars])")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:coefs] == $(m.internalModel.bounding_constr_mip[i][:coefs])")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:sense] == :($(m.internalModel.bounding_constr_mip[i][:sense]))")
    #         println("@test m.internalModel.bounding_constr_mip[$(i)][:cnt] == $(m.internalModel.bounding_constr_mip[i][:cnt])")
    #     else
    #         error("Unknown structural_constr type $(m.structural_constr[i])")
    #     end
    # end

    @testset "Expression Hanlding Corner Cases - 1 : sign convertor special case" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log_level=100)

        m = Model(solver=test_solver)
        @variable(m, x[1:5], Bin)
        @NLconstraint(m, x[1] + -x[2] >= 2)
        @objective(m, Min, x[1]+x[2])
        JuMP.build(m)

        @test m.internalModel.bounding_constr_mip[1][:rhs] == 2.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, -1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 2
    end

    @testset "Expression Hanlding Corner Cases - 2 : full sub-expression" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log_level=100)

        m = Model(solver=test_solver)
        @variable(m, x[1:5]>=0)
        @NLconstraint(m, (x[1] + 10 + 4 * x[2] + 200) * 5 >= 2 * 5)
        @constraint(m, x[1] + x[2] + 200 >= 2)
        @objective(m, Min, x[1]+x[2])
        JuMP.build(m)

        @test m.internalModel.bounding_constr_mip[1][:rhs] == -198.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[2][:rhs] == -1040.0
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[5.0, 20.0]
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 2
    end

end
