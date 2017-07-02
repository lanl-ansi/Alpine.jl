@testset "Expression Parser Test" begin

    @testset "Expression Test || bilinear || Affine || exprs.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0)

        m=exprstest(solver=test_solver)

        JuMP.build(m)

        ex = m.internalModel.lifted_constr_expr_mip[1]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [-1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[1][:coefs]
        @test affdict[:vars] == [:(x[1])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[1][:vars]
        @test isapprox(affdict[:rhs], 109.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[1][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[1][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[2]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [3.0,3.0,3.0,3.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[2][:coefs]
        @test affdict[:vars] == [:(x[8]),:(x[9]),:(x[10]),:(x[11])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[2][:vars]
        @test isapprox(affdict[:rhs], 111.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[2][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[2][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[3]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [-1.0,20.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[3][:coefs]
        @test affdict[:vars] == [:(x[12]),:(x[13])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[3][:vars]
        @test isapprox(affdict[:rhs], 222.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[3][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[3][:sense]

        # 1.0 * x[12] - 115.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[4]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [-1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[4][:coefs]
        @test affdict[:vars] == [:(x[12])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[4][:vars]
        @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[4][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[4][:sense]

        # 1.0 * x[12] - 115.0 <= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[5]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[5][:coefs]
        @test affdict[:vars] == [:(x[12])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[5][:vars]
        @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[5][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[5][:sense]

        # -1.0 * x[12] - 115.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[6]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [-1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[6][:coefs]
        @test affdict[:vars] == [:(x[12])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[6][:vars]
        @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[6][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[6][:sense]

        # (x[1] + 1.0 * x[14]) - 555.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[7]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [1.0, 1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[7][:coefs]
        @test affdict[:vars] == [:(x[1]),:(x[14])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[7][:vars]
        @test isapprox(affdict[:rhs], 555.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[7][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[7][:sense]

        # ((x[8] - 7.0 * x[9]) + x[10] + x[4]) - 6666.0 <= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[8]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [1.0,-7.0,1.0,1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[8][:coefs]
        @test affdict[:vars] == [:(x[8]),:(x[9]),:(x[10]),:(x[4])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[8][:vars]
        @test isapprox(affdict[:rhs], 6666.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[8][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[8][:sense]

        # ((13.0 * x[1] - x[2]) + 30.0 * x[3] + x[4]) - 77.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[9]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [13.0,-1.0,30.0,1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[9][:coefs]
        @test affdict[:vars] == [:(x[1]),:(x[2]),:(x[3]),:(x[4])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[9][:vars]
        @test isapprox(affdict[:rhs], 77.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[9][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[9][:sense]
    end

    @testset "Expression Test || bilinear || Affine || nlp1.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0)
        m=nlp1(solver=test_solver)

        JuMP.build(m)

        ex = m.internalModel.lifted_constr_expr_mip[1]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [1.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[1][:coefs]
        @test affdict[:vars] == [:(x[5])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[1][:vars]
        @test isapprox(affdict[:rhs], 8.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[1][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[1][:sense]

    end

    @testset "Expression Test || bilinear || Affine || nlp3.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=CbcSolver(OutputFlag=0),
								   log_level=0)

        m=nlp3(solver=test_solver)

        JuMP.build(m)

        ex = m.internalModel.lifted_constr_expr_mip[1]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [0.0025,0.0025]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[1][:coefs]
        @test affdict[:vars] == [:(x[4]),:(x[6])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[1][:vars]
        @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[1][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[1][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[2]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [0.0025,-0.0025,0.0025]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[2][:coefs]
        @test affdict[:vars] == [:(x[5]),:(x[4]),:(x[7])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[2][:vars]
        @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[2][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[2][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[3]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [0.01, -0.01]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[3][:coefs]
        @test affdict[:vars] == [:(x[8]),:(x[5])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[3][:vars]
        @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[3][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[3][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[4]
        affdict = POD.expr_to_affine(ex)
        @test (affdict[:coefs] .== [100.0, -1.0, 833.33252]) == [true, true, true]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[4][:coefs]
        @test affdict[:vars] == [:(x[1]),:(x[9]),:(x[4])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[4][:vars]
        @test isapprox(affdict[:rhs], 83333.333; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[4][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[4][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[5]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [1.0,-1.0,-1250.0,1250.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[5][:coefs]
        @test affdict[:vars] == [:(x[10]),:(x[11]),:(x[4]),:(x[5])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[5][:vars]
        @test isapprox(affdict[:rhs], 0.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[5][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[5][:sense]

        ex = m.internalModel.lifted_constr_expr_mip[6]
        affdict = POD.expr_to_affine(ex)
        @test affdict[:coefs] == [1.0,-1.0,-2500.0]
        @test affdict[:coefs] == m.internalModel.lifted_constr_aff_mip[6][:coefs]
        @test affdict[:vars] == [:(x[12]),:(x[13]),:(x[5])]
        @test affdict[:vars] == m.internalModel.lifted_constr_aff_mip[6][:vars]
        @test isapprox(affdict[:rhs], -1.25e6; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.lifted_constr_aff_mip[6][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.lifted_constr_aff_mip[6][:sense]

	end

    @testset "Expression Test || bilinear || Simple || bi1.jl " begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0)

        m = operator_c(solver=test_solver)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonlinear_info)) == 8
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 1), Expr(:ref, :x, 1)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 2), Expr(:ref, :x, 2)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 3), Expr(:ref, :x, 3)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 4), Expr(:ref, :x, 4)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 2), Expr(:ref, :x, 3)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 3), Expr(:ref, :x, 4)])

        # TODO setup detailed check on this problem

    end
end
