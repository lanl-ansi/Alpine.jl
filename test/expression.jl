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

    @testset "Expression Test || bilinear || Complex || blend029.jl " begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0)

        m = blend029(solver=test_solver)

        JuMP.build(m) # Setup internal model
        @test length(keys(m.internalModel.nonlinear_info)) == 28
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 37), Expr(:ref, :x, 55)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 38), Expr(:ref, :x, 56)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 37), Expr(:ref, :x, 26)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 43), Expr(:ref, :x, 26)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 47), Expr(:ref, :x, 59)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 48), Expr(:ref, :x, 60)])
        @test haskey(m.internalModel.nonlinear_info, [Expr(:ref, :x, 47), Expr(:ref, :x, 36)])

        @test m.internalModel.lifted_constr_aff_mip[1][:rhs] == 1.0
        @test m.internalModel.lifted_constr_aff_mip[1][:vars] == Any[:(x[1]), :(x[4]), :(x[7]), :(x[10]), :(x[49])]
        @test m.internalModel.lifted_constr_aff_mip[1][:sense] == :(==)
        @test m.internalModel.lifted_constr_aff_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0, 1.0]
        @test m.internalModel.lifted_constr_aff_mip[1][:cnt] == 5

        @test m.internalModel.lifted_constr_aff_mip[4][:rhs] == 0.1
        @test m.internalModel.lifted_constr_aff_mip[4][:vars] == Any[:(x[4]), :(x[16]), :(x[25]), :(x[34]), :(x[58])]
        @test m.internalModel.lifted_constr_aff_mip[4][:sense] == :(==)
        @test m.internalModel.lifted_constr_aff_mip[4][:coefs] == Any[-1.0, -1.0, -1.0, 1.0, 1.0]
        @test m.internalModel.lifted_constr_aff_mip[4][:cnt] == 5

        @test m.internalModel.lifted_constr_aff_mip[17][:rhs] == -0.14
        @test m.internalModel.lifted_constr_aff_mip[17][:vars] == Any[:(x[11]), :(x[23]), :(x[32]), :(x[64]), :(x[65])]
        @test m.internalModel.lifted_constr_aff_mip[17][:sense] == :(==)
        @test m.internalModel.lifted_constr_aff_mip[17][:coefs] == Any[-1.0, -1.0, -1.0, -1.0, 1.0]
        @test m.internalModel.lifted_constr_aff_mip[17][:cnt] == 5

        @test m.internalModel.lifted_constr_aff_mip[67][:rhs] == 0.0
        @test m.internalModel.lifted_constr_aff_mip[67][:vars] == Any[:(x[13])]
        @test m.internalModel.lifted_constr_aff_mip[67][:sense] == :(>=)
        @test m.internalModel.lifted_constr_aff_mip[67][:coefs] == Any[1.0]
        @test m.internalModel.lifted_constr_aff_mip[67][:cnt] == 1

        @test m.internalModel.lifted_constr_aff_mip[103][:rhs] == 1.0
        @test m.internalModel.lifted_constr_aff_mip[103][:vars] == Any[:(x[73])]
        @test m.internalModel.lifted_constr_aff_mip[103][:sense] == :(<=)
        @test m.internalModel.lifted_constr_aff_mip[103][:coefs] == Any[1.0]
        @test m.internalModel.lifted_constr_aff_mip[103][:cnt] == 1

        @test m.internalModel.lifted_constr_aff_mip[127][:rhs] == -1.0
        @test m.internalModel.lifted_constr_aff_mip[127][:vars] == Any[:(x[73])]
        @test m.internalModel.lifted_constr_aff_mip[127][:sense] == :(>=)
        @test m.internalModel.lifted_constr_aff_mip[127][:coefs] == Any[-1.0]
        @test m.internalModel.lifted_constr_aff_mip[127][:cnt] == 1

        @test m.internalModel.lifted_constr_aff_mip[187][:rhs] == 1.0
        @test m.internalModel.lifted_constr_aff_mip[187][:vars] == Any[:(x[79]), :(x[94])]
        @test m.internalModel.lifted_constr_aff_mip[187][:sense] == :(<=)
        @test m.internalModel.lifted_constr_aff_mip[187][:coefs] == Any[1.0, 1.0]
        @test m.internalModel.lifted_constr_aff_mip[187][:cnt] == 2

        @test m.internalModel.lifted_constr_aff_mip[202][:rhs] == 0.04
        @test m.internalModel.lifted_constr_aff_mip[202][:vars] == Any[:(x[103]), :(x[1]), :(x[13]), :(x[25]), :(x[28]), :(x[31])]
        @test m.internalModel.lifted_constr_aff_mip[202][:sense] == :(==)
        @test m.internalModel.lifted_constr_aff_mip[202][:coefs] == Any[1.0, -0.6, -0.2, 0.2, 0.2, 0.2]
        @test m.internalModel.lifted_constr_aff_mip[202][:cnt] == 6

        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[55])]][:id] == 1
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[55])]][:lifted_var_ref] == :(x[103])
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[55])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[55])]][:nonlinear_type] == :bilinear

        @test m.internalModel.lifted_constr_aff_mip[206][:rhs] == 0.0
        @test m.internalModel.lifted_constr_aff_mip[206][:vars] == Any[:(x[107]), :(x[103]), :(x[108]), :(x[109]), :(x[110]), :(x[2]), :(x[14])]
        @test m.internalModel.lifted_constr_aff_mip[206][:sense] == :(==)
        @test m.internalModel.lifted_constr_aff_mip[206][:coefs] == Any[1.0, -1.0, 1.0, 1.0, 1.0, -0.6, -0.2]
        @test m.internalModel.lifted_constr_aff_mip[206][:cnt] == 7

        @test m.internalModel.nonlinear_info[[:(x[38]), :(x[56])]][:id] == 5
        @test m.internalModel.nonlinear_info[[:(x[38]), :(x[56])]][:lifted_var_ref] == :(x[107])
        @test m.internalModel.nonlinear_info[[:(x[38]), :(x[56])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[38]), :(x[56])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[26])]][:id] == 6
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[26])]][:lifted_var_ref] == :(x[108])
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[26])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[26])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[32])]][:id] == 8
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[32])]][:lifted_var_ref] == :(x[110])
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[32])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[37]), :(x[32])]][:nonlinear_type] == :bilinear

        @test m.internalModel.lifted_constr_aff_mip[213][:rhs] == 0.0
        @test m.internalModel.lifted_constr_aff_mip[213][:vars] == Any[:(x[129]), :(x[127]), :(x[124]), :(x[130]), :(x[6]), :(x[18])]
        @test m.internalModel.lifted_constr_aff_mip[213][:sense] == :(==)
        @test m.internalModel.lifted_constr_aff_mip[213][:coefs] == Any[1.0, -1.0, -1.0, 1.0, -0.4, -0.4]
        @test m.internalModel.lifted_constr_aff_mip[213][:cnt] == 6

        @test m.internalModel.nonlinear_info[[:(x[48]), :(x[60])]][:id] == 27
        @test m.internalModel.nonlinear_info[[:(x[48]), :(x[60])]][:lifted_var_ref] == :(x[129])
        @test m.internalModel.nonlinear_info[[:(x[48]), :(x[60])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[48]), :(x[60])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[59])]][:id] == 25
        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[59])]][:lifted_var_ref] == :(x[127])
        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[59])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[59])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[36])]][:id] == 28
        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[36])]][:lifted_var_ref] == :(x[130])
        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[36])]][:convexified] == false
        @test m.internalModel.nonlinear_info[[:(x[47]), :(x[36])]][:nonlinear_type] == :bilinear
    end

    @testset "Expression Test || multilinear || Simple || multi.jl " begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log_level=0)

        m = multi3(solver=test_solver, exprmode=1)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonlinear_info)) == 1
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[3]), :(x[2])])][:id] == 1
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[3]), :(x[2])])][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[3]), :(x[2])])][:nonlinear_type] == :multilinear

        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[4])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        @test m.internalModel.lifted_constr_aff_mip[1][:rhs] == 3.0
        @test m.internalModel.lifted_constr_aff_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3])]
        @test m.internalModel.lifted_constr_aff_mip[1][:sense] == :(<=)
        @test m.internalModel.lifted_constr_aff_mip[1][:coefs] == Any[1.0, 1.0, 1.0]
        @test m.internalModel.lifted_constr_aff_mip[1][:cnt] == 3

        m = multi3(solver=test_solver, exprmode=2)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear

        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[5])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi3(solver=test_solver, exprmode=3)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear

        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[5])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=1)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 1

        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[5])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[3]), :(x[2]), :(x[4])])][:id] == 1
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[3]), :(x[2]), :(x[4])])][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[3]), :(x[2]), :(x[4])])][:nonlinear_type] == :multilinear

        @test m.internalModel.lifted_constr_aff_mip[1][:rhs] == 4.0
        @test m.internalModel.lifted_constr_aff_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        @test m.internalModel.lifted_constr_aff_mip[1][:sense] == :(<=)
        @test m.internalModel.lifted_constr_aff_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0]
        @test m.internalModel.lifted_constr_aff_mip[1][:cnt] == 4

        m = multi4(solver=test_solver, exprmode=2)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 3

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:id] == 3
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:nonlinear_type] == :bilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing


        m = multi4(solver=test_solver, exprmode=3)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[Set(Any[:(x[5]), :(x[3]), :(x[4])])][:id] == 2
        @test m.internalModel.nonlinear_info[Set(Any[:(x[5]), :(x[3]), :(x[4])])][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[5]), :(x[3]), :(x[4])])][:nonlinear_type] == :multilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=4)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[Set(Any[:(x[5]), :(x[1]), :(x[2])])][:id] == 2
        @test m.internalModel.nonlinear_info[Set(Any[:(x[5]), :(x[1]), :(x[2])])][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[5]), :(x[1]), :(x[2])])][:nonlinear_type] == :multilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=5)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 3

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 3)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=6)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2

        @test m.internalModel.nonlinear_info[Set(Any[:(x[3]), :(x[1]), :(x[2])])][:id] == 1
        @test m.internalModel.nonlinear_info[Set(Any[:(x[3]), :(x[1]), :(x[2])])][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[3]), :(x[1]), :(x[2])])][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear

        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=7)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 3

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] == :bilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing


        m = multi4(solver=test_solver, exprmode=8)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[5]), :(x[4])])][:id] == 2
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[5]), :(x[4])])][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[1]), :(x[5]), :(x[4])])][:nonlinear_type] == :multilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing


        m = multi4(solver=test_solver, exprmode=9)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 2

        @test m.internalModel.nonlinear_info[Set(Any[:(x[3]), :(x[4]), :(x[2])])][:id] == 1
        @test m.internalModel.nonlinear_info[Set(Any[:(x[3]), :(x[4]), :(x[2])])][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[Set(Any[:(x[3]), :(x[4]), :(x[2])])][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear

        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=10)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 3

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] == :bilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=11)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_info)) == 3

        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_info[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear


        @test m.internalModel.lifted_obj_aff_mip[:rhs] == 0
        @test m.internalModel.lifted_obj_aff_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.lifted_obj_aff_mip[:coefs] == [1.0]
        @test m.internalModel.lifted_obj_aff_mip[:cnt] == 1
        @test m.internalModel.lifted_obj_aff_mip[:sense] == nothing

    end

end
