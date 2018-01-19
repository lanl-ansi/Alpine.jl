@testset "Expression Parser Test :: Basic I" begin

    @testset "Expression Test || bilinear || Affine || exprs.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log=10000)

        m=exprstest(solver=test_solver)

        JuMP.build(m)

        ex = m.internalModel.bounding_constr_expr_mip[1]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [-1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[1][:coefs]
        @test affdict[:vars] == [:(x[1])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[1][:vars]
        @test isapprox(affdict[:rhs], 109.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[1][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[1][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[2]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [3.0,3.0,3.0,3.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[2][:coefs]
        @test affdict[:vars] == [:(x[8]),:(x[9]),:(x[10]),:(x[11])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[2][:vars]
        @test isapprox(affdict[:rhs], 111.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[2][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[2][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[3]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [-1.0,20.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[3][:coefs]
        @test affdict[:vars] == [:(x[12]),:(x[13])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[3][:vars]
        @test isapprox(affdict[:rhs], 222.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[3][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[3][:sense]

        # 1.0 * x[12] - 115.0 >= 0.0
        ex = m.internalModel.bounding_constr_expr_mip[4]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [-1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[4][:coefs]
        @test affdict[:vars] == [:(x[12])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[4][:vars]
        @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[4][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[4][:sense]

        # 1.0 * x[12] - 115.0 <= 0.0
        ex = m.internalModel.bounding_constr_expr_mip[5]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[5][:coefs]
        @test affdict[:vars] == [:(x[12])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[5][:vars]
        @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[5][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[5][:sense]

        # -1.0 * x[12] - 115.0 >= 0.0
        ex = m.internalModel.bounding_constr_expr_mip[6]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [-1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[6][:coefs]
        @test affdict[:vars] == [:(x[12])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[6][:vars]
        @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[6][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[6][:sense]

        # (x[1] + 1.0 * x[14]) - 555.0 >= 0.0
        ex = m.internalModel.bounding_constr_expr_mip[7]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [1.0, 1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[7][:coefs]
        @test affdict[:vars] == [:(x[1]),:(x[14])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[7][:vars]
        @test isapprox(affdict[:rhs], 555.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[7][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[7][:sense]

        # ((x[8] - 7.0 * x[9]) + x[10] + x[4]) - 6666.0 <= 0.0
        ex = m.internalModel.bounding_constr_expr_mip[8]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [1.0,-7.0,1.0,1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[8][:coefs]
        @test affdict[:vars] == [:(x[8]),:(x[9]),:(x[10]),:(x[4])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[8][:vars]
        @test isapprox(affdict[:rhs], 6666.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[8][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[8][:sense]

        # ((13.0 * x[1] - x[2]) + 30.0 * x[3] + x[4]) - 77.0 >= 0.0
        ex = m.internalModel.bounding_constr_expr_mip[9]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [13.0,-1.0,30.0,1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[9][:coefs]
        @test affdict[:vars] == [:(x[1]),:(x[2]),:(x[3]),:(x[4])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[9][:vars]
        @test isapprox(affdict[:rhs], 77.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[9][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[9][:sense]
    end

    @testset "Expression Test || bilinear || Affine || nlp1.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000)
        m=nlp1(solver=test_solver)

        JuMP.build(m)

        ex = m.internalModel.bounding_constr_expr_mip[1]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [1.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[1][:coefs]
        @test affdict[:vars] == [:(x[5])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[1][:vars]
        @test isapprox(affdict[:rhs], 8.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[1][:rhs]
        @test affdict[:sense] == :(>=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[1][:sense]
    end

    @testset "Expression Test || bilinear || Affine || nlp3.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=CbcSolver(OutputFlag=0),
								   log=10000)

        m=nlp3(solver=test_solver)

        JuMP.build(m)

        ex = m.internalModel.bounding_constr_expr_mip[1]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [0.0025,0.0025]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[1][:coefs]
        @test affdict[:vars] == [:(x[4]),:(x[6])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[1][:vars]
        @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[1][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[1][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[2]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [0.0025,-0.0025,0.0025]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[2][:coefs]
        @test affdict[:vars] == [:(x[5]),:(x[4]),:(x[7])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[2][:vars]
        @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[2][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[2][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[3]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [0.01, -0.01]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[3][:coefs]
        @test affdict[:vars] == [:(x[8]),:(x[5])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[3][:vars]
        @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[3][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[3][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[4]
        affdict = POD.expr_linear_to_affine(ex)
        @test (affdict[:coefs] .== [100.0, -1.0, 833.33252]) == [true, true, true]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[4][:coefs]
        @test affdict[:vars] == [:(x[1]),:(x[9]),:(x[4])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[4][:vars]
        @test isapprox(affdict[:rhs], 83333.333; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[4][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[4][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[5]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [1.0,-1.0,-1250.0,1250.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[5][:coefs]
        @test affdict[:vars] == [:(x[10]),:(x[11]),:(x[4]),:(x[5])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[5][:vars]
        @test isapprox(affdict[:rhs], 0.0; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[5][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[5][:sense]

        ex = m.internalModel.bounding_constr_expr_mip[6]
        affdict = POD.expr_linear_to_affine(ex)
        @test affdict[:coefs] == [1.0,-1.0,-2500.0]
        @test affdict[:coefs] == m.internalModel.bounding_constr_mip[6][:coefs]
        @test affdict[:vars] == [:(x[12]),:(x[13]),:(x[5])]
        @test affdict[:vars] == m.internalModel.bounding_constr_mip[6][:vars]
        @test isapprox(affdict[:rhs], -1.25e6; atol = 1e-3)
        @test affdict[:rhs] == m.internalModel.bounding_constr_mip[6][:rhs]
        @test affdict[:sense] == :(<=)
        @test affdict[:sense] == m.internalModel.bounding_constr_mip[6][:sense]
	end

    @testset "Expression Test || bilinear || Simple || bi1.jl " begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000)

        m = operator_c(solver=test_solver)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonlinear_terms)) == 8
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 1), Expr(:ref, :x, 1)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 2)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 3)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 4), Expr(:ref, :x, 4)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 3)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 4)])

        # TODO setup detailed check on this problem
    end

    @testset "Expression Test || bilinear || Complex || blend029.jl " begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000)

        m = blend029(solver=test_solver)

        JuMP.build(m) # Setup internal model
        @test length(keys(m.internalModel.nonlinear_terms)) == 28
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 37), Expr(:ref, :x, 55)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 38), Expr(:ref, :x, 56)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 37), Expr(:ref, :x, 26)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 43), Expr(:ref, :x, 26)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 47), Expr(:ref, :x, 59)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 48), Expr(:ref, :x, 60)])
        @test haskey(m.internalModel.nonlinear_terms, [Expr(:ref, :x, 47), Expr(:ref, :x, 36)])

        @test m.internalModel.bounding_constr_mip[1][:rhs] == 1.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[4]), :(x[7]), :(x[10]), :(x[49])]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 5

        @test m.internalModel.bounding_constr_mip[4][:rhs] == 0.1
        @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[4]), :(x[16]), :(x[25]), :(x[34]), :(x[58])]
        @test m.internalModel.bounding_constr_mip[4][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[4][:coefs] == Any[-1.0, -1.0, -1.0, 1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[4][:cnt] == 5

        @test m.internalModel.bounding_constr_mip[17][:rhs] == -0.14
        @test m.internalModel.bounding_constr_mip[17][:vars] == Any[:(x[11]), :(x[23]), :(x[32]), :(x[64]), :(x[65])]
        @test m.internalModel.bounding_constr_mip[17][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[17][:coefs] == Any[-1.0, -1.0, -1.0, -1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[17][:cnt] == 5

        @test m.internalModel.bounding_constr_mip[67][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[67][:vars] == Any[:(x[13])]
        @test m.internalModel.bounding_constr_mip[67][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[67][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[67][:cnt] == 1

        @test m.internalModel.bounding_constr_mip[103][:rhs] == 1.0
        @test m.internalModel.bounding_constr_mip[103][:vars] == Any[:(x[73])]
        @test m.internalModel.bounding_constr_mip[103][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[103][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[103][:cnt] == 1

        @test m.internalModel.bounding_constr_mip[127][:rhs] == -1.0
        @test m.internalModel.bounding_constr_mip[127][:vars] == Any[:(x[73])]
        @test m.internalModel.bounding_constr_mip[127][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[127][:coefs] == Any[-1.0]
        @test m.internalModel.bounding_constr_mip[127][:cnt] == 1

        @test m.internalModel.bounding_constr_mip[187][:rhs] == 1.0
        @test m.internalModel.bounding_constr_mip[187][:vars] == Any[:(x[79]), :(x[94])]
        @test m.internalModel.bounding_constr_mip[187][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[187][:coefs] == Any[1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[187][:cnt] == 2

        @test m.internalModel.bounding_constr_mip[202][:rhs] == 0.04
        @test m.internalModel.bounding_constr_mip[202][:vars] == Any[:(x[103]), :(x[1]), :(x[13]), :(x[25]), :(x[28]), :(x[31])]
        @test m.internalModel.bounding_constr_mip[202][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[202][:coefs] == Any[1.0, -0.6, -0.2, 0.2, 0.2, 0.2]
        @test m.internalModel.bounding_constr_mip[202][:cnt] == 6

        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[55])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[55])]][:lifted_var_ref] == :(x[103])
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[55])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[55])]][:nonlinear_type] == :bilinear

        @test m.internalModel.bounding_constr_mip[206][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[206][:vars] == Any[:(x[107]), :(x[103]), :(x[108]), :(x[109]), :(x[110]), :(x[2]), :(x[14])]
        @test m.internalModel.bounding_constr_mip[206][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[206][:coefs] == Any[1.0, -1.0, 1.0, 1.0, 1.0, -0.6, -0.2]
        @test m.internalModel.bounding_constr_mip[206][:cnt] == 7

        @test m.internalModel.nonlinear_terms[[:(x[38]), :(x[56])]][:id] == 5
        @test m.internalModel.nonlinear_terms[[:(x[38]), :(x[56])]][:lifted_var_ref] == :(x[107])
        @test m.internalModel.nonlinear_terms[[:(x[38]), :(x[56])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[38]), :(x[56])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[26])]][:id] == 6
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[26])]][:lifted_var_ref] == :(x[108])
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[26])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[26])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[32])]][:id] == 8
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[32])]][:lifted_var_ref] == :(x[110])
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[32])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[37]), :(x[32])]][:nonlinear_type] == :bilinear

        @test m.internalModel.bounding_constr_mip[213][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[213][:vars] == Any[:(x[129]), :(x[127]), :(x[124]), :(x[130]), :(x[6]), :(x[18])]
        @test m.internalModel.bounding_constr_mip[213][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[213][:coefs] == Any[1.0, -1.0, -1.0, 1.0, -0.4, -0.4]
        @test m.internalModel.bounding_constr_mip[213][:cnt] == 6

        @test m.internalModel.nonlinear_terms[[:(x[48]), :(x[60])]][:id] == 27
        @test m.internalModel.nonlinear_terms[[:(x[48]), :(x[60])]][:lifted_var_ref] == :(x[129])
        @test m.internalModel.nonlinear_terms[[:(x[48]), :(x[60])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[48]), :(x[60])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[59])]][:id] == 25
        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[59])]][:lifted_var_ref] == :(x[127])
        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[59])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[59])]][:nonlinear_type] == :bilinear

        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[36])]][:id] == 28
        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[36])]][:lifted_var_ref] == :(x[130])
        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[36])]][:convexified] == false
        @test m.internalModel.nonlinear_terms[[:(x[47]), :(x[36])]][:nonlinear_type] == :bilinear
    end

    @testset "Expression Test || multilinear || Simple || multi.jl " begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000)

        m = multi3(solver=test_solver, exprmode=1)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonlinear_terms)) == 1
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3])]][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :multilinear

        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[4])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        @test m.internalModel.bounding_constr_mip[1][:rhs] == 3.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3])]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 3

        m = multi3(solver=test_solver, exprmode=2)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear

        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[5])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi3(solver=test_solver, exprmode=3)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear

        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[5])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=1)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 1

        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[5])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :multilinear

        @test m.internalModel.bounding_constr_mip[1][:rhs] == 4.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 4

        m = multi4(solver=test_solver, exprmode=2)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 3

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:id] == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:nonlinear_type] == :bilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing


        m = multi4(solver=test_solver, exprmode=3)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[3]), :(x[4])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[3]), :(x[4])]][:nonlinear_type] == :multilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=4)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[5])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[5])]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[5])]][:nonlinear_type] == :multilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=5)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 3

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=6)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3])]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear

        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=7)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 3

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] == :bilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing


        m = multi4(solver=test_solver, exprmode=8)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[5]), :(x[4])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[5]), :(x[4])]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[5]), :(x[4])]][:nonlinear_type] == :multilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing


        m = multi4(solver=test_solver, exprmode=9)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 2

        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[3]), :(x[4])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear

        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=10)

        JuMP.build(m)
        @test length(keys(m.internalModel.nonlinear_terms)) == 3

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] == :bilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing

        m = multi4(solver=test_solver, exprmode=11)

        JuMP.build(m)

        @test length(keys(m.internalModel.nonlinear_terms)) == 3

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[7])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] == :bilinear


        @test m.internalModel.bounding_obj_mip[:rhs] == 0
        @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
        @test m.internalModel.bounding_obj_mip[:cnt] == 1
        @test m.internalModel.bounding_obj_mip[:sense] == nothing
    end

    @testset "Expression Test || bilinear || Complex-div || div.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000)

        m = div(solver=test_solver)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonlinear_terms)) == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:id] == 1
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:lifted_var_ref] == :(x[3])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:nonlinear_type] == :monomial

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:id] == 2
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[4])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:nonlinear_type] == :monomial

        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 3
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
        @test m.internalModel.nonlinear_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :bilinear

        aff_mip = m.internalModel.bounding_constr_mip

        @test aff_mip[1][:rhs] == 0.0
        @test aff_mip[1][:vars] == Any[:(x[1])]
        @test aff_mip[1][:sense] == :(>=)
        @test round.(aff_mip[1][:coefs],1) == Any[-3.0]
        @test aff_mip[1][:cnt] == 1

        @test aff_mip[2][:rhs] == 0.0
        @test aff_mip[2][:vars] == Any[:(x[2])]
        @test aff_mip[2][:sense] == :(>=)
        @test aff_mip[2][:coefs] == Any[0.25]
        @test aff_mip[2][:cnt] == 1

        @test aff_mip[3][:rhs] == 0.0
        @test aff_mip[3][:vars] == Any[:(x[2])]
        @test aff_mip[3][:sense] == :(>=)
        @test aff_mip[3][:coefs] == Any[5.0]
        @test aff_mip[3][:cnt] == 1

        @test aff_mip[4][:rhs] == 0.0
        @test aff_mip[4][:vars] == Any[:(x[2])]
        @test aff_mip[4][:sense] == :(>=)
        @test aff_mip[4][:coefs] == Any[-120.0]
        @test aff_mip[4][:cnt] == 1

        @test aff_mip[5][:rhs] == 0.0
        @test aff_mip[5][:vars] == Any[:(x[2])]
        @test aff_mip[5][:sense] == :(>=)
        @test aff_mip[5][:coefs] == Any[72000.0]
        @test aff_mip[5][:cnt] == 1

        @test aff_mip[6][:rhs] == 0.0
        @test aff_mip[6][:vars] == Any[:(x[1])]
        @test aff_mip[6][:sense] == :(>=)
        @test aff_mip[6][:coefs] == Any[72000.0]
        @test aff_mip[6][:cnt] == 1

        @test aff_mip[7][:rhs] == 8.0
        @test aff_mip[7][:vars] == Any[:(x[5])]
        @test aff_mip[7][:sense] == :(>=)
        @test aff_mip[7][:coefs] == Any[0.6]
        @test aff_mip[7][:cnt] == 1

        @test aff_mip[8][:rhs] == 0.0
        @test aff_mip[8][:vars] == Any[:(x[2]),:(x[5])]
        @test aff_mip[8][:sense] == :(>=)
        @test aff_mip[8][:coefs] == Any[5.6,-72000.0]
        @test aff_mip[8][:cnt] == 2

        @test aff_mip[9][:rhs] == 0.0
        @test aff_mip[9][:vars] == Any[:(x[2]),:(x[5])]
        @test aff_mip[9][:sense] == :(>=)
        @test aff_mip[9][:coefs] == Any[5.6,-36000.0]
        @test aff_mip[9][:cnt] == 2

        @test aff_mip[10][:rhs] == 0.0
        @test aff_mip[10][:vars] == Any[:(x[2]),:(x[5]),:(x[2])]
        @test aff_mip[10][:sense] == :(>=)
        @test aff_mip[10][:coefs] == Any[5.6, -300.0, -1.75]
        @test aff_mip[10][:cnt] == 3

        @test aff_mip[11][:rhs] == 0.0
        @test aff_mip[11][:vars] == Any[:(x[2]),:(x[1]),:(x[5])]
        @test aff_mip[11][:sense] == :(>=)
        @test aff_mip[11][:coefs] == Any[5.6, -0.5, 0.5]
        @test aff_mip[11][:cnt] == 3
    end
end

@testset "Expression Parser Test :: Basic II" begin

    @testset "Expression parsing || part1 " begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(OutputFlag=0),
								   log=10000))
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

    @testset "Expression parsing || part2" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000))

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
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[5])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[1])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[2])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]), :(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[9]), :(x[10])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[8]), :(x[10])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[12])]) #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[7])])  #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[14]), :(x[10])]) #15
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[3])]) #16
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]), :(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[8])]) #18
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[18])]) #19

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[5])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:id] == 3
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:id] == 4
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[8])]][:id] == 5
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[3])]][:id] == 6
        @test m.internalModel.nonlinear_terms[[:(x[9]), :(x[10])]][:id] == 7
        @test m.internalModel.nonlinear_terms[[:(x[8]), :(x[10])]][:id] == 8
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[12])]][:id] == 9
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[7])]][:id] == 10
        @test m.internalModel.nonlinear_terms[[:(x[14]), :(x[10])]][:id] == 11
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[3])]][:id] == 12
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[16])]][:id] == 13
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[8])]][:id] == 14
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[18])]][:id] == 15
    end

    @testset "Expression parsing || part3" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000))

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
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[5])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[6])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[1])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[9])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[10])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[2])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[12])])  #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[13])]) #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[14])]) #15
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]), :(x[3])]) #16
        @test haskey(m.internalModel.nonlinear_terms, [:(x[5]), :(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[17])]) #18
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[4])]) #19
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]), :(x[19])]) #20
        @test haskey(m.internalModel.nonlinear_terms, [:(x[8]), :(x[12])]) #21
        @test haskey(m.internalModel.nonlinear_terms, [:(x[21]), :(x[16])]) #22
        @test haskey(m.internalModel.nonlinear_terms, [:(x[22]), :(x[19])]) #23

        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[2])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[5])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[6])]][:id] == 3
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:id] == 4
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[8])]][:id] == 5
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[9])]][:id] == 6
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[10])]][:id] == 7
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:id] == 8
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[12])]][:id] == 9
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[13])]][:id] == 10
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[14])]][:id] == 11
        @test m.internalModel.nonlinear_terms[[:(x[3]), :(x[3])]][:id] == 12
        @test m.internalModel.nonlinear_terms[[:(x[5]), :(x[16])]][:id] == 13
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[17])]][:id] == 14
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[4])]][:id] == 15
        @test m.internalModel.nonlinear_terms[[:(x[6]), :(x[19])]][:id] == 16
        @test m.internalModel.nonlinear_terms[[:(x[8]), :(x[12])]][:id] == 17
        @test m.internalModel.nonlinear_terms[[:(x[21]), :(x[16])]][:id] == 18
        @test m.internalModel.nonlinear_terms[[:(x[22]), :(x[19])]][:id] == 19
    end

    @testset "Expression parsing || part7" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000))
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]*x[2]*x[3]^2*x[4]^2 >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]^2*x[4] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]*x[4]^2 >= 1)
        @NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4]^2 >= 1)

        JuMP.build(m)

        @test m.internalModel.bounding_constr_expr_mip[1] == :(x[5]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[2] == :(x[7]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[3] == :(x[9]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[4] == :(x[12]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[5] == :(x[13]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[6] == :(x[14]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[7] == :(x[15]-1.0 >= 0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2]),:(x[3]),:(x[4])]) #5
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[1])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[3]),:(x[4]),:(x[6])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[2])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[3]),:(x[4]),:(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[3])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[4])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2]),:(x[10]),:(x[11])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[4]),:(x[8]),:(x[10])]) #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[3]),:(x[6]),:(x[11])]) #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[6]),:(x[8]),:(x[10]),:(x[11])]) #15

        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[2]),:(x[3]),:(x[4])]][:id] == 1 #5
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[1])]][:id] == 2 #6
        @test m.internalModel.nonlinear_terms[[:(x[2]),:(x[3]),:(x[4]),:(x[6])]][:id] == 3 #7
        @test m.internalModel.nonlinear_terms[[:(x[2]),:(x[2])]][:id] == 4 #8
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[3]),:(x[4]),:(x[8])]][:id] == 5 #9
        @test m.internalModel.nonlinear_terms[[:(x[3]),:(x[3])]][:id] == 6 #10
        @test m.internalModel.nonlinear_terms[[:(x[4]),:(x[4])]][:id] == 7 #11
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[2]),:(x[10]),:(x[11])]][:id] ==  8  #12
        @test m.internalModel.nonlinear_terms[[:(x[1]),:(x[4]),:(x[8]),:(x[10])]][:id] == 9 #13
        @test m.internalModel.nonlinear_terms[[:(x[2]),:(x[3]),:(x[6]),:(x[11])]][:id] == 10 #14
        @test m.internalModel.nonlinear_terms[[:(x[6]),:(x[8]),:(x[10]),:(x[11])]][:id] == 11 #15
    end

    @testset "Expression parsing || part8" begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000))
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

        @test m.internalModel.bounding_constr_expr_mip[1] == :(x[6]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[2] == :(x[9]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[3] == :(x[12]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[4] == :(x[15]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[5] == :(x[17]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[6] == :(x[19]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[7] == :(x[21]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[8] == :(x[22]-1.0 >= 0.0)
        @test m.internalModel.bounding_constr_expr_mip[9] == :(x[24]-1.0 >= 0.0)

        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2]),:(x[3])]) #5
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[5])]) #6
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[1])]) #7
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[3]),:(x[7])]) #8
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[8])]) #9
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[2])]) #10
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[10])]) #11
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[4]),:(x[11])]) #12
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[3])]) #13
        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]),:(x[13])]) #14
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[4]),:(x[14])]) #15
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[10])]) #16
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[4]),:(x[16])]) #17
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]),:(x[2])]) #18
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[18]),:(x[13])]) #19
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]),:(x[4])]) #20
        @test haskey(m.internalModel.nonlinear_terms, [:(x[3]),:(x[18]),:(x[20])]) #21
        @test haskey(m.internalModel.nonlinear_terms, [:(x[18]),:(x[13]),:(x[20])]) #22
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]),:(x[10]),:(x[13])]) #23
        @test haskey(m.internalModel.nonlinear_terms, [:(x[23]),:(x[20])]) #24
    end
end

@testset "Expression Parsing Test :: Convex" begin

    @testset "Convex Parsing :: PART I" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log=10000)
        m = convex_test(test_solver)

        JuMP.build(m)

        @test m.internalModel.num_constr_convex == 21

        # 0 : OBJ
        @test m.internalModel.structural_obj == :convex
        @test m.internalModel.nonlinear_constrs[0][:expr_orig] == :objective
        @test m.internalModel.nonlinear_constrs[0][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[0][:convexified] == false

        @test m.internalModel.bounding_obj_mip[:sense] == nothing
        @test m.internalModel.bounding_obj_mip[:coefs] == [1.0, 1.0]
        @test m.internalModel.bounding_obj_mip[:vars] == [:(x[1]), :(x[3])]
        @test m.internalModel.bounding_obj_mip[:rhs] == 0.0
        @test m.internalModel.bounding_obj_mip[:powers] == [2, 2]
        @test m.internalModel.bounding_obj_mip[:cnt] == 2

        # 1
        @test m.internalModel.structural_constr[1] == :convex
        @test m.internalModel.nonlinear_constrs[1][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[1][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[1][:convexified] == false
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[1][:coefs] == [3.0, 4.0]
        @test m.internalModel.bounding_constr_mip[1][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[1][:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 2

        # 2
        @test m.internalModel.structural_constr[2] == :convex
        @test m.internalModel.nonlinear_constrs[2][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[2][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[2][:convexified] == false
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[2][:coefs] == [3.0, 4.0]
        @test m.internalModel.bounding_constr_mip[2][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[2][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[2][:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 2

        # 4
        @test m.internalModel.structural_constr[4] == :convex
        @test m.internalModel.nonlinear_constrs[4][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[4][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[4][:convexified] == false
        @test m.internalModel.bounding_constr_mip[4][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[4][:coefs] == [3.0, 4.0]
        @test m.internalModel.bounding_constr_mip[4][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[4][:rhs] == 10.0
        @test m.internalModel.bounding_constr_mip[4][:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[4][:cnt] == 2

        # 5
        @test m.internalModel.structural_constr[5] == :convex
        @test m.internalModel.nonlinear_constrs[5][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[5][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[5][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[5][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[5][:coefs] == [3.0, 4.0, 6.0]
        @test m.internalModel.bounding_constr_mip[5][:vars] == [:(x[1]), :(x[2]), :(x[3])]
        @test m.internalModel.bounding_constr_mip[5][:rhs] == 10.0
        @test m.internalModel.bounding_constr_mip[5][:powers] == [2, 2, 2]
        @test m.internalModel.bounding_constr_mip[5][:cnt] == 3

        # 6
        @test m.internalModel.structural_constr[6] == :convex
        @test m.internalModel.nonlinear_constrs[6][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[6][:convex_type] == :convexC
        @test m.internalModel.nonlinear_constrs[6][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[6][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[6][:coefs] == [3.0, 4.0, 5.0]
        @test m.internalModel.bounding_constr_mip[6][:vars] == [:(x[1]), :(x[2]), :(x[5])]
        @test m.internalModel.bounding_constr_mip[6][:rhs] == 100.0
        @test m.internalModel.bounding_constr_mip[6][:powers] == [0.5, 0.5, 0.5]
        @test m.internalModel.bounding_constr_mip[6][:cnt] == 3

        # 7
        @test m.internalModel.structural_constr[7] == :convex
        @test m.internalModel.nonlinear_constrs[7][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[7][:convex_type] == :convexC
        @test m.internalModel.nonlinear_constrs[7][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[7][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[7][:coefs] == [-3.0, -4.0]
        @test m.internalModel.bounding_constr_mip[7][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[7][:rhs] == -100.0
        @test m.internalModel.bounding_constr_mip[7][:powers] == [0.5, 0.5]
        @test m.internalModel.bounding_constr_mip[7][:cnt] == 2

        # 8
        @test m.internalModel.structural_constr[8] == :convex
        @test m.internalModel.nonlinear_constrs[8][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[8][:convex_type] == :convexB
        @test m.internalModel.nonlinear_constrs[8][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[8][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[8][:coefs] == [3.0, 1.0, 5.0]
        @test m.internalModel.bounding_constr_mip[8][:vars] == [:(x[1]), :(x[2]), :(x[3])]
        @test m.internalModel.bounding_constr_mip[8][:rhs] == 200.0
        @test m.internalModel.bounding_constr_mip[8][:powers] == [3.0, 3.0, 3.0]
        @test m.internalModel.bounding_constr_mip[8][:cnt] == 3

        # 9
        @test m.internalModel.structural_constr[9] == :convex
        @test m.internalModel.nonlinear_constrs[9][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[9][:convex_type] == :convexB
        @test m.internalModel.nonlinear_constrs[9][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[9][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[9][:coefs] == [1.0, 1.0, 1.0, 100.0]
        @test m.internalModel.bounding_constr_mip[9][:vars] == [:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        @test m.internalModel.bounding_constr_mip[9][:rhs] == 200.0
        @test m.internalModel.bounding_constr_mip[9][:powers] == [3, 3, 3, 3]
        @test m.internalModel.bounding_constr_mip[9][:cnt] == 4

        # 11
        @test m.internalModel.structural_constr[11] == :convex
        @test m.internalModel.nonlinear_constrs[11][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[11][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[11][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[11][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[11][:coefs] == [3.0, 4.0]
        @test m.internalModel.bounding_constr_mip[11][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[11][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[11][:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[11][:cnt] == 2

        # 14
        @test m.internalModel.structural_constr[14] == :convex
        @test m.internalModel.nonlinear_constrs[14][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[14][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[14][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[14][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[14][:coefs] == [3.0, 5.0]
        @test m.internalModel.bounding_constr_mip[14][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[14][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[14][:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[14][:cnt] == 2

        # 15
        @test m.internalModel.structural_constr[15] == :convex
        @test m.internalModel.nonlinear_constrs[15][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[15][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[15][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[15][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[15][:coefs] == [3.0, 5.0, 1.0]
        @test m.internalModel.bounding_constr_mip[15][:vars] == [:(x[1]), :(x[2]), :(x[4])]
        @test m.internalModel.bounding_constr_mip[15][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[15][:powers] == [2, 2, 2]
        @test m.internalModel.bounding_constr_mip[15][:cnt] == 3

        # 19
        @test m.internalModel.structural_constr[19] == :convex
        @test m.internalModel.nonlinear_constrs[19][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[19][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[19][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[19][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[19][:coefs] == [3.0, 16.0]
        @test m.internalModel.bounding_constr_mip[19][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[19][:rhs] == 40.0
        @test m.internalModel.bounding_constr_mip[19][:powers] == [2, 2]
        @test m.internalModel.bounding_constr_mip[19][:cnt] == 2

        # 22
        @test m.internalModel.structural_constr[22] == :convex
        @test m.internalModel.nonlinear_constrs[22][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[22][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[22][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[22][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[22][:coefs] == [3.0, 4.0, 5.0, 6.0]
        @test m.internalModel.bounding_constr_mip[22][:vars] == [:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        @test m.internalModel.bounding_constr_mip[22][:rhs] == 15.0
        @test m.internalModel.bounding_constr_mip[22][:powers] == [2, 2, 2, 2]
        @test m.internalModel.bounding_constr_mip[22][:cnt] == 4

        # 25
        @test m.internalModel.structural_constr[25] == :convex
        @test m.internalModel.nonlinear_constrs[25][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[25][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[25][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[25][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[25][:coefs] == [1.0, 1.0, 1.0, 1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[25][:vars] == [:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]
        @test m.internalModel.bounding_constr_mip[25][:rhs] == 99999.0
        @test m.internalModel.bounding_constr_mip[25][:powers] == [2, 2, 2, 2, 2]
        @test m.internalModel.bounding_constr_mip[25][:cnt] == 5

        # 26
        @test m.internalModel.structural_constr[26] == :convex
        @test m.internalModel.nonlinear_constrs[26][:expr_orig] == :constraints
        @test m.internalModel.nonlinear_constrs[26][:convex_type] == :convexA
        @test m.internalModel.nonlinear_constrs[26][:convexified] == :false
        @test m.internalModel.bounding_constr_mip[26][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[26][:coefs] == [3.0, 4.0]
        @test m.internalModel.bounding_constr_mip[26][:vars] == [:(x[1]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[26][:rhs] == 200.0
        @test m.internalModel.bounding_constr_mip[26][:powers] == [4, 4]
        @test m.internalModel.bounding_constr_mip[26][:cnt] == 2
    end
end

@testset "Expression Praser Test :: Linear Lifting" begin

    @testset "Linear Lifting : nlp2" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=PajaritoSolver(mip_solver=CbcSolver(logLevel=0),cont_solver=IpoptSolver(print_level=0), log_level=0),
                               disc_ratio=8,
                               log=10000)

        m = nlp2(solver=test_solver)

        JuMP.build(m)

        @test length(m.internalModel.linear_terms) == 2
        @test length(m.internalModel.nonlinear_terms) == 4

        lk = Vector{Any}(2)
        lk[1] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, -1.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3)])))
        lk[2] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, -2.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6)])))

        @test length(keys(m.internalModel.linear_terms)) == 2
        yids = [4,7]
        for i in 1:length(keys(m.internalModel.linear_terms))
            for j in keys(m.internalModel.linear_terms)
                if m.internalModel.linear_terms[j][:id] == i
                    @test j == lk[i]
                    @test j[:sign] == :+
                    @test j[:scalar] == lk[i][:scalar]
                    @test j[:coef_var] == lk[i][:coef_var]
                    @test m.internalModel.linear_terms[j][:y_idx] == yids[i]
                end
            end
        end

        @test haskey(m.internalModel.linear_terms, lk[1])
        @test haskey(m.internalModel.linear_terms, lk[2])


        @test haskey(m.internalModel.nonlinear_terms, [:(x[2]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[4]), :(x[4])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[7]), :(x[7])])
        @test haskey(m.internalModel.nonlinear_terms, [:(x[1]), :(x[1])])
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:id] == 1
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:id] == 3
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[4])]][:id] == 2
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[7])]][:id] == 4
        @test m.internalModel.nonlinear_terms[[:(x[1]), :(x[1])]][:lifted_var_ref].args[2] == 3
        @test m.internalModel.nonlinear_terms[[:(x[2]), :(x[2])]][:lifted_var_ref].args[2] == 6
        @test m.internalModel.nonlinear_terms[[:(x[4]), :(x[4])]][:lifted_var_ref].args[2] == 5
        @test m.internalModel.nonlinear_terms[[:(x[7]), :(x[7])]][:lifted_var_ref].args[2] == 8
    end

    @testset "Linear Lifting : general" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(OutputFlag=0),
                               log=10000)

        m = basic_linear_lift(solver=test_solver)

        JuMP.build(m) # Setup internal model

        lk = Vector{Any}(5)
        lk[1] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1), (-1.0, 2)])))
        lk[2] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(3.0, 2), (-1.0, 3)])))
        lk[3] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1), (1.0, 3)])))
        lk[4] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1)])))
        lk[5] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 3.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1), (1.0, 3)])))

        @test length(keys(m.internalModel.linear_terms)) == 5
        yids = [8,9,12,14,16]
        for i in 1:length(keys(m.internalModel.linear_terms))
            for j in keys(m.internalModel.linear_terms)
                if m.internalModel.linear_terms[j][:id] == i
                    @test j == lk[i]
                    @test j[:sign] == :+
                    @test j[:scalar] == lk[i][:scalar]
                    @test j[:coef_var] == lk[i][:coef_var]
                    @test m.internalModel.linear_terms[j][:y_idx] == yids[i]
                end
            end
        end

        nlk1 = [:(x[8]), :(x[9]), :(x[12])]
        nlk2 = [:(x[2]), :(x[2])]
        nlk3 = [:(x[2]), :(x[3])]
        nlk4 = [:(x[8]), :(x[8])]
        nlk5 = [:(x[1]), :(x[3])]
        nlk6 = [:(x[8]), :(x[9])]
        nlk7 = [:(x[1]), :(x[2])]
        nlk8 = [:(x[16]), :(x[15])]
        nlk9 = [:(x[14]), :(x[14])]

        @test m.internalModel.nonlinear_terms[nlk1][:id] == 7
        @test m.internalModel.nonlinear_terms[nlk2][:id] == 3
        @test m.internalModel.nonlinear_terms[nlk3][:id] == 4
        @test m.internalModel.nonlinear_terms[nlk4][:id] == 6
        @test m.internalModel.nonlinear_terms[nlk5][:id] == 2
        @test m.internalModel.nonlinear_terms[nlk6][:id] == 5
        @test m.internalModel.nonlinear_terms[nlk7][:id] == 1
        @test m.internalModel.nonlinear_terms[nlk8][:id] == 9
        @test m.internalModel.nonlinear_terms[nlk9][:id] == 8

        @test m.internalModel.nonlinear_terms[nlk1][:lifted_var_ref].args[2] == 13
        @test m.internalModel.nonlinear_terms[nlk2][:lifted_var_ref].args[2] == 6
        @test m.internalModel.nonlinear_terms[nlk3][:lifted_var_ref].args[2] == 7
        @test m.internalModel.nonlinear_terms[nlk4][:lifted_var_ref].args[2] == 11
        @test m.internalModel.nonlinear_terms[nlk5][:lifted_var_ref].args[2] == 5
        @test m.internalModel.nonlinear_terms[nlk6][:lifted_var_ref].args[2] == 10
        @test m.internalModel.nonlinear_terms[nlk7][:lifted_var_ref].args[2] == 4
        @test m.internalModel.nonlinear_terms[nlk8][:lifted_var_ref].args[2] == 17
        @test m.internalModel.nonlinear_terms[nlk9][:lifted_var_ref].args[2] == 15

        @test m.internalModel.nonlinear_terms[nlk1][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_terms[nlk2][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk3][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk4][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk5][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk6][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk7][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk8][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk9][:nonlinear_type] == :monomial

        @test m.internalModel.nonlinear_terms[nlk1][:var_idxs] == [8,9,12]
        @test m.internalModel.nonlinear_terms[nlk2][:var_idxs] == [2,2]
        @test m.internalModel.nonlinear_terms[nlk3][:var_idxs] == [2,3]
        @test m.internalModel.nonlinear_terms[nlk4][:var_idxs] == [8]
        @test m.internalModel.nonlinear_terms[nlk5][:var_idxs] == [1,3]
        @test m.internalModel.nonlinear_terms[nlk6][:var_idxs] == [8,9]
        @test m.internalModel.nonlinear_terms[nlk7][:var_idxs] == [1,2]
        @test m.internalModel.nonlinear_terms[nlk8][:var_idxs] == [16, 15]
        @test m.internalModel.nonlinear_terms[nlk9][:var_idxs] == [14]
    end

    @testset "Expression Test || complex || Affine || operator_b" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log=10000)

        m=operator_b(solver=test_solver)

        JuMP.build(m)

        nlk = Vector{Any}(31)
        nlk[1] = [:(x[1]), :(x[2])]
        nlk[2] = [:(x[2]), :(x[3])]
        nlk[3] = [:(x[1]), :(x[10])]
        nlk[4] = [:(x[2]), :(x[12])]
        nlk[5] = [:(x[1]), :(x[1])]
        nlk[6] = [:(x[2]), :(x[2])]
        nlk[7] = [:(x[3]), :(x[3]), :(x[3])]
        nlk[8] = [:(x[17]), :(x[17])]
        nlk[9] = [:(x[1]), :(x[1]), :(x[1])]
        nlk[10] = Expr[:(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3]), :(x[3])]
        nlk[11] = [:(x[1]), :(x[2]), :(x[3])]
        nlk[12] = [:(x[3]), :(x[8])]
        nlk[13] = [:(x[1]), :(x[9])]
        nlk[14] = [:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        nlk[15] = [:(x[3]), :(x[4])]
        nlk[16] = [:(x[8]), :(x[25])]
        nlk[17] = [:(x[1]), :(x[4]), :(x[9])]
        nlk[18] = [:(x[3]), :(x[4]), :(x[8])]
        nlk[19] = [:(x[4]), :(x[22])]
        nlk[20] = [:(x[3]), :(x[3])]
        nlk[21] = [:(x[4]), :(x[14]), :(x[15]), :(x[30])]
        nlk[22] = [:(x[1]), :(x[1]), :(x[2]), :(x[2]), :(x[3]), :(x[3])]
        nlk[23] = [:(x[33]), :(x[34]), :(x[35]), :(x[36])]
        nlk[24] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
        nlk[25] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
        nlk[26] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[40]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
        nlk[27] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[42]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
        nlk[28] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[44]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
        nlk[29] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[46]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
        nlk[30] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[48]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
        nlk[31] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[21]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))

        @test length(keys(m.internalModel.nonlinear_terms)) == 31
        ids = [1:31;]
        yids = [8,9,11,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,37,38,39,41,43,45,47,49,50]
        nltypes = [:bilinear, :bilinear, :bilinear, :bilinear, :monomial,
                   :monomial, :multilinear, :monomial, :multilinear, :multilinear,
                   :multilinear, :bilinear, :bilinear, :multilinear, :bilinear,
                   :bilinear, :multilinear, :multilinear, :bilinear, :monomial,
                   :multilinear, :multilinear, :multilinear, :sin, :cos,
                   :sin, :cos, :sin, :sin, :cos, :sin]
        for i in 1:31
            @test m.internalModel.nonlinear_terms[nlk[i]][:id] == i
            @test m.internalModel.nonlinear_terms[nlk[i]][:y_idx] == yids[i]
            @test m.internalModel.nonlinear_terms[nlk[i]][:nonlinear_type] == nltypes[i]
            @test m.internalModel.nonlinear_terms[nlk[i]][:y_type] == :Cont
        end

        lk = Vector{Any}(12)
        lk[1] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 4.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2)])))
        lk[2] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 3)])))
        lk[3] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 5.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1)])))
        lk[4] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1)])))
        lk[5] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 2.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2)])))
        lk[6] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 3.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3)])))
        lk[7] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 4.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4)])))
        lk[8] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1)])))
        lk[9] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (-1.0, 3)])))
        lk[10] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 2), (1.0, 3)])))
        lk[11] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(4.0, 8)])))
        lk[12] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 9)])))

        @test length(keys(m.internalModel.linear_terms)) == 12
        yids = [10,12,17,33,34,35,36,40,42,44,46,48]
        scas = [5.0,1.0,6.0,2.0,3.0,4.0,5.0,1.0,1.0,1.0,0.0,0.0]
        for i in 1:length(keys(m.internalModel.linear_terms))
            for j in keys(m.internalModel.linear_terms)
                if m.internalModel.linear_terms[j][:id] == i
                    @test j == lk[i]
                    @test j[:sign] == :+
                    @test j[:scalar] == lk[i][:scalar]
                    @test j[:coef_var] == lk[i][:coef_var]
                    @test m.internalModel.linear_terms[j][:y_idx] == yids[i]
                end
            end
        end
    end
end

@testset "Expression Parser Test :: BMPL -> BINLIN Opt" begin

    @testset "Operator :: bmpl && linbin && binprod" begin

        test_solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(logLevel=0),
                               log=10000)

        m = bpml(solver=test_solver)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonlinear_terms)) == 12

        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[6])]][:y_idx] == 11
        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 12
        @test m.internalModel.nonlinear_terms[Expr[:(x[12]), :(x[6])]][:y_idx] == 13
        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:y_idx] == 14
        @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[6])]][:y_idx] == 15
        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7])]][:y_idx] == 16
        @test m.internalModel.nonlinear_terms[Expr[:(x[16]), :(x[1])]][:y_idx] == 17
        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7]), :(x[8])]][:y_idx] == 18
        @test m.internalModel.nonlinear_terms[Expr[:(x[18]), :(x[1])]][:y_idx] == 19
        @test m.internalModel.nonlinear_terms[Expr[:(x[18]), :(x[14])]][:y_idx] == 20
        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]][:y_idx] == 21
        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7]), :(x[9]), :(x[10]), :(x[6])]][:y_idx] == 22

        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[Expr[:(x[12]), :(x[6])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[6])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[16]), :(x[1])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7]), :(x[8])]][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[18]), :(x[1])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[18]), :(x[14])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7]), :(x[9]), :(x[10]), :(x[6])]][:nonlinear_type] == :multilinear

        @test length(m.internalModel.var_type_lifted) == 22
        @test m.internalModel.var_type_lifted[11] == :Cont
        @test m.internalModel.var_type_lifted[12] == :Bin
        @test m.internalModel.var_type_lifted[13] == :Cont
        @test m.internalModel.var_type_lifted[14] == :Bin
        @test m.internalModel.var_type_lifted[15] == :Cont
        @test m.internalModel.var_type_lifted[16] == :Cont
        @test m.internalModel.var_type_lifted[17] == :Cont
        @test m.internalModel.var_type_lifted[18] == :Cont
        @test m.internalModel.var_type_lifted[19] == :Cont
        @test m.internalModel.var_type_lifted[20] == :Cont
        @test m.internalModel.var_type_lifted[21] == :Bin
        @test m.internalModel.var_type_lifted[22] == :Cont
    end

    @testset "Operator :: bmpl && linbin && binprod with linear lifting and coefficients" begin
        test_solver=PODSolver(nlp_local_solver=IpoptSolver(),
                               mip_solver=CbcSolver(logLevel=0),
                               log=10000)

        m = bmpl_linearlifting(solver=test_solver)

        JuMP.build(m)

        lk1 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),
                               Pair{Symbol,Any}(:scalar, 0.0),
                               Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1), (1.0, 17)])))
        lk2 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),
                               Pair{Symbol,Any}(:scalar, 0.0),
                               Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 11), (1.0, 12)])))
        lk3 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),
                               Pair{Symbol,Any}(:scalar, 0.0),
                               Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7), (2.0, 6)])))
        lk4 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),
                               Pair{Symbol,Any}(:scalar, 0.0),
                               Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6), (1.0, 19)])))
        lk5 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),
                               Pair{Symbol,Any}(:scalar, 0.0),
                               Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6),  (1.0, 25)])))

        @test haskey(m.internalModel.linear_terms, lk1)
        @test haskey(m.internalModel.linear_terms, lk2)
        @test haskey(m.internalModel.linear_terms, lk3)
        @test haskey(m.internalModel.linear_terms, lk4)
        @test haskey(m.internalModel.linear_terms, lk5)

        @test m.internalModel.linear_terms[lk1][:lifted_constr_ref] == :(x[20] == x[1] + x[17])
        @test m.internalModel.linear_terms[lk2][:lifted_constr_ref] == :(x[13] == x[11] + x[12])
        @test m.internalModel.linear_terms[lk3][:lifted_constr_ref] == :(x[15] == 2*x[6] + x[7])
        @test m.internalModel.linear_terms[lk4][:lifted_constr_ref] == :(x[21] == x[6] + x[19])
        @test m.internalModel.linear_terms[lk5][:lifted_constr_ref] == :(x[26] == x[6] + x[25])

        @test m.internalModel.linear_terms[lk1][:y_idx] == 20
        @test m.internalModel.linear_terms[lk2][:y_idx] == 13
        @test m.internalModel.linear_terms[lk3][:y_idx] == 15
        @test m.internalModel.linear_terms[lk4][:y_idx] == 21
        @test m.internalModel.linear_terms[lk5][:y_idx] == 26

        @test m.internalModel.linear_terms[lk1][:y_type] == :Cont
        @test m.internalModel.linear_terms[lk2][:y_type] == :Cont
        @test m.internalModel.linear_terms[lk3][:y_type] == :Cont
        @test m.internalModel.linear_terms[lk4][:y_type] == :Cont
        @test m.internalModel.linear_terms[lk5][:y_type] == :Cont

        nlk1 = Expr[:(x[18]), :(x[7])]
        nlk2 = Expr[:(x[6]), :(x[7]), :(x[8])]
        nlk3 = Expr[:(x[7]), :(x[10])]
        nlk4 = Expr[:(x[2]), :(x[26])]
        nlk5 = Expr[:(x[1]), :(x[13])]
        nlk6 = Expr[:(x[1]), :(x[15])]
        nlk7 = Expr[:(x[2]), :(x[6])]
        nlk8 = Expr[:(x[23]), :(x[24])]
        nlk9 = Expr[:(x[1]), :(x[2])]
        nlk10 = Expr[:(x[2]), :(x[3])]
        nlk11 = Expr[:(x[20]), :(x[21])]
        nlk12 = Expr[:(x[3]), :(x[4])]

        @test haskey(m.internalModel.nonlinear_terms, nlk1)
        @test haskey(m.internalModel.nonlinear_terms, nlk2)
        @test haskey(m.internalModel.nonlinear_terms, nlk3)
        @test haskey(m.internalModel.nonlinear_terms, nlk4)
        @test haskey(m.internalModel.nonlinear_terms, nlk5)
        @test haskey(m.internalModel.nonlinear_terms, nlk6)
        @test haskey(m.internalModel.nonlinear_terms, nlk7)
        @test haskey(m.internalModel.nonlinear_terms, nlk8)
        @test haskey(m.internalModel.nonlinear_terms, nlk9)
        @test haskey(m.internalModel.nonlinear_terms, nlk10)
        @test haskey(m.internalModel.nonlinear_terms, nlk11)
        @test haskey(m.internalModel.nonlinear_terms, nlk12)


        @test m.internalModel.nonlinear_terms[nlk1][:id] == 7
        @test m.internalModel.nonlinear_terms[nlk2][:id] == 1
        @test m.internalModel.nonlinear_terms[nlk3][:id] == 9
        @test m.internalModel.nonlinear_terms[nlk4][:id] == 12
        @test m.internalModel.nonlinear_terms[nlk5][:id] == 3
        @test m.internalModel.nonlinear_terms[nlk6][:id] == 4
        @test m.internalModel.nonlinear_terms[nlk7][:id] == 5
        @test m.internalModel.nonlinear_terms[nlk8][:id] == 11
        @test m.internalModel.nonlinear_terms[nlk9][:id] == 6
        @test m.internalModel.nonlinear_terms[nlk10][:id] == 2
        @test m.internalModel.nonlinear_terms[nlk11][:id] == 8
        @test m.internalModel.nonlinear_terms[nlk12][:id] == 10


        @test m.internalModel.nonlinear_terms[nlk1][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk2][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk3][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk4][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk5][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk6][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk7][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk8][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk9][:y_type] == :Bin
        @test m.internalModel.nonlinear_terms[nlk10][:y_type] == :Bin
        @test m.internalModel.nonlinear_terms[nlk11][:y_type] == :Cont
        @test m.internalModel.nonlinear_terms[nlk12][:y_type] == :Bin

        @test m.internalModel.nonlinear_terms[nlk1][:y_idx] == 19
        @test m.internalModel.nonlinear_terms[nlk2][:y_idx] == 11
        @test m.internalModel.nonlinear_terms[nlk3][:y_idx] == 23
        @test m.internalModel.nonlinear_terms[nlk4][:y_idx] == 27
        @test m.internalModel.nonlinear_terms[nlk5][:y_idx] == 14
        @test m.internalModel.nonlinear_terms[nlk6][:y_idx] == 16
        @test m.internalModel.nonlinear_terms[nlk7][:y_idx] == 17
        @test m.internalModel.nonlinear_terms[nlk8][:y_idx] == 25
        @test m.internalModel.nonlinear_terms[nlk9][:y_idx] == 18
        @test m.internalModel.nonlinear_terms[nlk10][:y_idx] == 12
        @test m.internalModel.nonlinear_terms[nlk11][:y_idx] == 22
        @test m.internalModel.nonlinear_terms[nlk12][:y_idx] == 24

        @test m.internalModel.nonlinear_terms[nlk1][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk2][:nonlinear_type] == :multilinear
        @test m.internalModel.nonlinear_terms[nlk3][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk4][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[nlk5][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[nlk6][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[nlk7][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk8][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk9][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[nlk10][:nonlinear_type] == :binprod
        @test m.internalModel.nonlinear_terms[nlk11][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[nlk12][:nonlinear_type] == :binprod
    end
end

@testset "Expression Parser Test :: Special Operators" begin
    test_solver=PODSolver(nlp_local_solver=IpoptSolver(),
                           mip_solver=CbcSolver(logLevel=0),
                           log=10000)

    m = specialopts(solver=test_solver)
    JuMP.build(m)

    nlk = Vector{Any}(10)
    nlk[1] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
    nlk[2] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
    nlk[3] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[13]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
    nlk[4] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[15]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
    nlk[5] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[11]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
    nlk[6] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[19]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))
    nlk[7] = Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[17]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))
    nlk[8] = Expr[:(x[7]), :(x[8])]
    nlk[9] = Expr[:(x[2]), :(x[7])]
    nlk[10] = Expr[:(x[1]), :(x[2])]

    for k in nlk
        @test haskey(m.internalModel.nonlinear_terms, k)
        isa(k, Dict) && @test k[:operator] == m.internalModel.nonlinear_terms[k][:nonlinear_type]
    end

    lk = Vector{Any}(4)
    lk[1] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.2, 2), (0.2, 1)])))
    lk[2] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.5, 1)])))
    lk[3] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.2, 11)])))
    lk[4] = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1)])))

    for k in lk
        @test haskey(m.internalModel.linear_terms, k)
    end

    @test m.internalModel.linear_terms[lk[1]][:lifted_constr_ref] == :(x[19] == 0.2 * x[1] + 0.2 * x[2])
    @test m.internalModel.linear_terms[lk[2]][:lifted_constr_ref] == :(x[15] == 0.5 * x[1])
    @test m.internalModel.linear_terms[lk[3]][:lifted_constr_ref] == :(x[17] == 0.2 * x[11])
    @test m.internalModel.linear_terms[lk[4]][:lifted_constr_ref] == :(x[13] == x[1] + x[2])
end

@testset "Expression Parser Test :: Basic Multiplication Operators (Machine Generated for diffs)" begin

    test_solver=PODSolver(nlp_local_solver=IpoptSolver(),
                           mip_solver=CbcSolver(logLevel=0),
                           log=10000)

    m = operator_basic(solver=test_solver)
    JuMP.build(m)

    # Test Generation Code commented
    # for k in keys(m.internalModel.nonlinear_terms)
    #     println("@test m.internalModel.nonlinear_terms[$(k)][:y_idx] == $(m.internalModel.nonlinear_terms[k][:y_idx])")
    #     println("@test m.internalModel.nonlinear_terms[$(k)][:id] == $(m.internalModel.nonlinear_terms[k][:id])")
    #     println("@test m.internalModel.nonlinear_terms[$(k)][:lifted_constr_ref] == :($(m.internalModel.nonlinear_terms[k][:lifted_constr_ref]))")
    #     println("@test m.internalModel.nonlinear_terms[$(k)][:nonlinear_type] == :$(m.internalModel.nonlinear_terms[k][:nonlinear_type])")
    # end

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

    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4])]][:y_idx] == 31
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4])]][:id] == 27
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[31] == x[3] * x[4])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:y_idx] == 83
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:id] == 79
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:lifted_constr_ref] == :(x[83] == x[4] * x[5] * x[76])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[2])]][:y_idx] == 9
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[2])]][:id] == 5
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[2])]][:lifted_constr_ref] == :(x[9] == (*)(x[2]))
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[2])]][:nonlinear_type] == :monomial
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7])]][:y_idx] == 19
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7])]][:id] == 15
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7])]][:lifted_constr_ref] == :(x[19] == x[5] * x[7])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[31])]][:y_idx] == 32
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[31])]][:id] == 28
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[31])]][:lifted_constr_ref] == :(x[32] == x[2] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[31])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:y_idx] == 49
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:id] == 45
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:lifted_constr_ref] == :(x[49] == x[1] * x[2] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:y_idx] == 50
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:id] == 46
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:lifted_constr_ref] == :(x[50] == x[3] * x[5] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[47])]][:y_idx] == 69
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[47])]][:id] == 65
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[47])]][:lifted_constr_ref] == :(x[69] == x[4] * x[47])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[47])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[4])]][:y_idx] == 28
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[4])]][:id] == 24
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[4])]][:lifted_constr_ref] == :(x[28] == (*)(x[4]))
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[4])]][:nonlinear_type] == :monomial
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[10])]][:y_idx] == 37
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[10])]][:id] == 33
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[10])]][:lifted_constr_ref] == :(x[37] == x[4] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[10])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[1])]][:y_idx] == 5
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[1])]][:id] == 1
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[1])]][:lifted_constr_ref] == :(x[5] == (*)(x[1]))
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[1])]][:nonlinear_type] == :monomial
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[95])]][:y_idx] == 96
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[95])]][:id] == 92
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[95])]][:lifted_constr_ref] == :(x[96] == x[5] * x[95])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[95])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[40])]][:y_idx] == 41
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[40])]][:id] == 37
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[40])]][:lifted_constr_ref] == :(x[41] == x[2] * x[40])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[40])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10])]][:y_idx] == 26
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10])]][:id] == 22
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10])]][:lifted_constr_ref] == :(x[26] == x[6] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:y_idx] == 84
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:id] == 80
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:lifted_constr_ref] == :(x[84] == x[1] * x[78] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[43])]][:y_idx] == 58
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[43])]][:id] == 54
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[43])]][:lifted_constr_ref] == :(x[58] == x[6] * x[43])
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[43])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:y_idx] == 106
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:id] == 102
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:lifted_constr_ref] == :(x[106] == x[4] * x[6] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:y_idx] == 90
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:id] == 86
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:lifted_constr_ref] == :(x[90] == x[2] * x[4] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[41])]][:y_idx] == 42
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[41])]][:id] == 38
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[41])]][:lifted_constr_ref] == :(x[42] == x[1] * x[41])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[41])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[37])]][:y_idx] == 38
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[37])]][:id] == 34
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[37])]][:lifted_constr_ref] == :(x[38] == x[2] * x[37])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[37])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[85])]][:y_idx] == 87
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[85])]][:id] == 83
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[85])]][:lifted_constr_ref] == :(x[87] == x[5] * x[85])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[85])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:y_idx] == 66
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:id] == 62
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:lifted_constr_ref] == :(x[66] == x[2] * x[3] * x[5] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:y_idx] == 53
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:id] == 49
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[53] == x[5] * x[9] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:y_idx] == 107
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:id] == 103
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:lifted_constr_ref] == :(x[107] == x[3] * x[6] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5])]][:y_idx] == 17
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5])]][:id] == 13
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5])]][:lifted_constr_ref] == :(x[17] == x[2] * x[5])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[37])]][:y_idx] == 57
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[37])]][:id] == 53
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[37])]][:lifted_constr_ref] == :(x[57] == x[6] * x[37])
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[37])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:y_idx] == 81
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:id] == 77
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:lifted_constr_ref] == :(x[81] == x[5] * x[7] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[38])]][:y_idx] == 39
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[38])]][:id] == 35
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[38])]][:lifted_constr_ref] == :(x[39] == x[1] * x[38])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[38])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:y_idx] == 48
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:id] == 44
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:lifted_constr_ref] == :(x[48] == x[1] * x[3] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[23])]][:y_idx] == 24
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[23])]][:id] == 20
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[23])]][:lifted_constr_ref] == :(x[24] == x[4] * x[23])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[23])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[17])]][:y_idx] == 23
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[17])]][:id] == 19
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[17])]][:lifted_constr_ref] == :(x[23] == x[3] * x[17])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[17])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:y_idx] == 51
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:id] == 47
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[51] == x[1] * x[9] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:y_idx] == 75
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:id] == 71
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:lifted_constr_ref] == :(x[75] == x[4] * x[5] * x[7])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:y_idx] == 95
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:id] == 91
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[95] == x[9] * x[10] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10])]][:y_idx] == 11
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10])]][:id] == 7
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[11] == x[9] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:y_idx] == 100
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:id] == 96
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:lifted_constr_ref] == :(x[100] == x[1] * x[2] * x[37])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:y_idx] == 62
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:id] == 58
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:lifted_constr_ref] == :(x[62] == x[2] * x[3] * x[4] * x[5])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:y_idx] == 65
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:id] == 61
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[65] == x[1] * x[4] * x[9] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:y_idx] == 80
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:id] == 76
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:lifted_constr_ref] == :(x[80] == x[1] * x[7] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:y_idx] == 46
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:id] == 42
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:lifted_constr_ref] == :(x[46] == x[1] * x[2] * x[3])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[40])]][:y_idx] == 59
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[40])]][:id] == 55
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[40])]][:lifted_constr_ref] == :(x[59] == x[17] * x[40])
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[40])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[85])]][:y_idx] == 86
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[85])]][:id] == 82
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[85])]][:lifted_constr_ref] == :(x[86] == x[1] * x[85])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[85])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:y_idx] == 97
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:id] == 93
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:lifted_constr_ref] == :(x[97] == x[1] * x[2] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[13]), :(x[28])]][:y_idx] == 29
    @test m.internalModel.nonlinear_terms[Expr[:(x[13]), :(x[28])]][:id] == 25
    @test m.internalModel.nonlinear_terms[Expr[:(x[13]), :(x[28])]][:lifted_constr_ref] == :(x[29] == x[13] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[13]), :(x[28])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:y_idx] == 92
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:id] == 88
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:lifted_constr_ref] == :(x[92] == x[2] * x[3] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9])]][:y_idx] == 20
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9])]][:id] == 16
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9])]][:lifted_constr_ref] == :(x[20] == x[1] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:y_idx] == 109
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:id] == 105
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[109] == x[14] * x[10] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[92])]][:y_idx] == 93
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[92])]][:id] == 89
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[92])]][:lifted_constr_ref] == :(x[93] == x[1] * x[92])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[92])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:y_idx] == 85
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:id] == 81
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[85] == x[2] * x[3] * x[4])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[11])]][:y_idx] == 16
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[11])]][:id] == 12
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[11])]][:lifted_constr_ref] == :(x[16] == x[1] * x[11])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[11])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[21])]][:y_idx] == 25
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[21])]][:id] == 21
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[21])]][:lifted_constr_ref] == :(x[25] == x[4] * x[21])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[21])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[53]), :(x[28])]][:y_idx] == 73
    @test m.internalModel.nonlinear_terms[Expr[:(x[53]), :(x[28])]][:id] == 69
    @test m.internalModel.nonlinear_terms[Expr[:(x[53]), :(x[28])]][:lifted_constr_ref] == :(x[73] == x[53] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[53]), :(x[28])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[32])]][:y_idx] == 34
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[32])]][:id] == 30
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[32])]][:lifted_constr_ref] == :(x[34] == x[5] * x[32])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[32])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[28])]][:y_idx] == 40
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[28])]][:id] == 36
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[28])]][:lifted_constr_ref] == :(x[40] == x[3] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[28])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:y_idx] == 67
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:id] == 63
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[67] == x[5] * x[9] * x[10] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:y_idx] == 102
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:id] == 98
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:lifted_constr_ref] == :(x[102] == x[5] * x[9] * x[43])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[26])]][:y_idx] == 27
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[26])]][:id] == 23
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[26])]][:lifted_constr_ref] == :(x[27] == x[4] * x[26])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[26])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[20]), :(x[31])]][:y_idx] == 56
    @test m.internalModel.nonlinear_terms[Expr[:(x[20]), :(x[31])]][:id] == 52
    @test m.internalModel.nonlinear_terms[Expr[:(x[20]), :(x[31])]][:lifted_constr_ref] == :(x[56] == x[20] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[20]), :(x[31])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:y_idx] == 63
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:id] == 59
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:lifted_constr_ref] == :(x[63] == x[1] * x[3] * x[4] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10])]][:y_idx] == 15
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10])]][:id] == 11
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10])]][:lifted_constr_ref] == :(x[15] == x[14] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[14]), :(x[10])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:y_idx] == 77
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:id] == 73
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:lifted_constr_ref] == :(x[77] == x[1] * x[4] * x[76])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[88])]][:y_idx] == 89
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[88])]][:id] == 85
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[88])]][:lifted_constr_ref] == :(x[89] == x[1] * x[88])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[88])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[90])]][:y_idx] == 94
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[90])]][:id] == 90
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[90])]][:lifted_constr_ref] == :(x[94] == x[5] * x[90])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[90])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:y_idx] == 108
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:id] == 104
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[108] == x[6] * x[10] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[48])]][:y_idx] == 70
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[48])]][:id] == 66
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[48])]][:lifted_constr_ref] == :(x[70] == x[4] * x[48])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[48])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[31])]][:y_idx] == 54
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[31])]][:id] == 50
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[31])]][:lifted_constr_ref] == :(x[54] == x[6] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[31])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:y_idx] == 99
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:id] == 95
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:lifted_constr_ref] == :(x[99] == x[1] * x[9] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[35])]][:y_idx] == 36
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[35])]][:id] == 32
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[35])]][:lifted_constr_ref] == :(x[36] == x[1] * x[35])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[35])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9])]][:y_idx] == 14
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9])]][:id] == 10
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9])]][:lifted_constr_ref] == :(x[14] == x[5] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[9])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[13])]][:y_idx] == 22
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[13])]][:id] == 18
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[13])]][:lifted_constr_ref] == :(x[22] == x[4] * x[13])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[13])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:y_idx] == 64
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:id] == 60
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[64] == x[1] * x[2] * x[10] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[46]), :(x[28])]][:y_idx] == 72
    @test m.internalModel.nonlinear_terms[Expr[:(x[46]), :(x[28])]][:id] == 68
    @test m.internalModel.nonlinear_terms[Expr[:(x[46]), :(x[28])]][:lifted_constr_ref] == :(x[72] == x[46] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[46]), :(x[28])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[31])]][:y_idx] == 55
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[31])]][:id] == 51
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[31])]][:lifted_constr_ref] == :(x[55] == x[17] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[31])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:y_idx] == 98
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:id] == 94
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:lifted_constr_ref] == :(x[98] == x[2] * x[5] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[10]), :(x[28])]][:y_idx] == 43
    @test m.internalModel.nonlinear_terms[Expr[:(x[10]), :(x[28])]][:id] == 39
    @test m.internalModel.nonlinear_terms[Expr[:(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[43] == x[10] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[10]), :(x[28])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[10])]][:y_idx] == 78
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[10])]][:id] == 74
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[10])]][:lifted_constr_ref] == :(x[78] == x[2] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[10])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:y_idx] == 52
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:id] == 48
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:lifted_constr_ref] == :(x[52] == x[2] * x[5] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[90])]][:y_idx] == 91
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[90])]][:id] == 87
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[90])]][:lifted_constr_ref] == :(x[91] == x[1] * x[90])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[90])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:y_idx] == 104
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:id] == 100
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:lifted_constr_ref] == :(x[104] == x[3] * x[4] * x[17])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[9])]][:y_idx] == 76
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[9])]][:id] == 72
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[9])]][:lifted_constr_ref] == :(x[76] == x[3] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[9])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:y_idx] == 79
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:id] == 75
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:lifted_constr_ref] == :(x[79] == x[1] * x[4] * x[78])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[3])]][:y_idx] == 10
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[3])]][:id] == 6
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[3])]][:lifted_constr_ref] == :(x[10] == (*)(x[3]))
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[3])]][:nonlinear_type] == :monomial
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[49])]][:y_idx] == 71
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[49])]][:id] == 67
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[49])]][:lifted_constr_ref] == :(x[71] == x[4] * x[49])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[49])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:y_idx] == 88
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:id] == 84
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:lifted_constr_ref] == :(x[88] == x[3] * x[4] * x[9])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[15]), :(x[28])]][:y_idx] == 30
    @test m.internalModel.nonlinear_terms[Expr[:(x[15]), :(x[28])]][:id] == 26
    @test m.internalModel.nonlinear_terms[Expr[:(x[15]), :(x[28])]][:lifted_constr_ref] == :(x[30] == x[15] * x[28])
    @test m.internalModel.nonlinear_terms[Expr[:(x[15]), :(x[28])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[46])]][:y_idx] == 68
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[46])]][:id] == 64
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[46])]][:lifted_constr_ref] == :(x[68] == x[4] * x[46])
    @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[46])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[11])]][:y_idx] == 12
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[11])]][:id] == 8
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[11])]][:lifted_constr_ref] == :(x[12] == x[5] * x[11])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[11])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:y_idx] == 74
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:id] == 70
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:lifted_constr_ref] == :(x[74] == x[1] * x[4] * x[7])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:y_idx] == 105
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:id] == 101
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:lifted_constr_ref] == :(x[105] == x[3] * x[4] * x[20])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:y_idx] == 82
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:id] == 78
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:lifted_constr_ref] == :(x[82] == x[1] * x[4] * x[11])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[43])]][:y_idx] == 44
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[43])]][:id] == 40
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[43])]][:lifted_constr_ref] == :(x[44] == x[9] * x[43])
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[43])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:y_idx] == 61
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 57
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[61] == x[1] * x[2] * x[3] * x[4])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[10])]][:y_idx] == 18
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[10])]][:id] == 14
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[10])]][:lifted_constr_ref] == :(x[18] == x[17] * x[10])
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[10])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7])]][:y_idx] == 8
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7])]][:id] == 4
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7])]][:lifted_constr_ref] == :(x[8] == x[1] * x[7])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[7])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 6
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2])]][:id] == 2
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2])]][:lifted_constr_ref] == :(x[6] == x[1] * x[2])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:y_idx] == 101
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:id] == 97
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:lifted_constr_ref] == :(x[101] == x[1] * x[2] * x[40])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[37])]][:y_idx] == 60
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[37])]][:id] == 56
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[37])]][:lifted_constr_ref] == :(x[60] == x[17] * x[37])
    @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[37])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:y_idx] == 103
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:id] == 99
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:lifted_constr_ref] == :(x[103] == x[3] * x[4] * x[6])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3])]][:y_idx] == 7
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3])]][:id] == 3
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3])]][:lifted_constr_ref] == :(x[7] == x[2] * x[3])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:y_idx] == 47
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:id] == 43
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:lifted_constr_ref] == :(x[47] == x[2] * x[3] * x[5])
    @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:nonlinear_type] == :multilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6])]][:y_idx] == 13
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6])]][:id] == 9
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6])]][:lifted_constr_ref] == :(x[13] == x[3] * x[6])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[6])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[31])]][:y_idx] == 35
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[31])]][:id] == 31
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[31])]][:lifted_constr_ref] == :(x[35] == x[9] * x[31])
    @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[31])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[32])]][:y_idx] == 33
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[32])]][:id] == 29
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[32])]][:lifted_constr_ref] == :(x[33] == x[1] * x[32])
    @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[32])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[20])]][:y_idx] == 21
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[20])]][:id] == 17
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[20])]][:lifted_constr_ref] == :(x[21] == x[3] * x[20])
    @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[20])]][:nonlinear_type] == :bilinear
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[44])]][:y_idx] == 45
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[44])]][:id] == 41
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[44])]][:lifted_constr_ref] == :(x[45] == x[5] * x[44])
    @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[44])]][:nonlinear_type] == :bilinear

    @test m.internalModel.bounding_constr_mip[1][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[5])]
    @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[1][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[1][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[2][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[6])]
    @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[2][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[3][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[5]), :(x[7])]
    @test m.internalModel.bounding_constr_mip[3][:coefs] == Any[1.0, 1.0]
    @test m.internalModel.bounding_constr_mip[3][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[3][:cnt] == 2
    @test m.internalModel.bounding_constr_mip[4][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[8])]
    @test m.internalModel.bounding_constr_mip[4][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[4][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[4][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[5][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[5][:vars] == Any[:(x[12])]
    @test m.internalModel.bounding_constr_mip[5][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[5][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[5][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[6][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[6][:vars] == Any[:(x[13])]
    @test m.internalModel.bounding_constr_mip[6][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[6][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[6][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[7][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[7][:vars] == Any[:(x[15])]
    @test m.internalModel.bounding_constr_mip[7][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[7][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[7][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[8][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[8][:vars] == Any[:(x[16])]
    @test m.internalModel.bounding_constr_mip[8][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[8][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[8][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[9][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[9][:vars] == Any[:(x[18])]
    @test m.internalModel.bounding_constr_mip[9][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[9][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[9][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[10][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[10][:vars] == Any[:(x[19])]
    @test m.internalModel.bounding_constr_mip[10][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[10][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[10][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[11][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[11][:vars] == Any[:(x[21])]
    @test m.internalModel.bounding_constr_mip[11][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[11][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[11][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[12][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[12][:vars] == Any[:(x[22])]
    @test m.internalModel.bounding_constr_mip[12][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[12][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[12][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[13][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[13][:vars] == Any[:(x[24])]
    @test m.internalModel.bounding_constr_mip[13][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[13][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[13][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[14][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[14][:vars] == Any[:(x[25])]
    @test m.internalModel.bounding_constr_mip[14][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[14][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[14][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[15][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[15][:vars] == Any[:(x[27])]
    @test m.internalModel.bounding_constr_mip[15][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[15][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[15][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[16][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[16][:vars] == Any[:(x[29])]
    @test m.internalModel.bounding_constr_mip[16][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[16][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[16][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[17][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[17][:vars] == Any[:(x[30])]
    @test m.internalModel.bounding_constr_mip[17][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[17][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[17][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[18][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[18][:vars] == Any[:(x[33])]
    @test m.internalModel.bounding_constr_mip[18][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[18][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[18][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[19][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[19][:vars] == Any[:(x[34])]
    @test m.internalModel.bounding_constr_mip[19][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[19][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[19][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[20][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[20][:vars] == Any[:(x[36])]
    @test m.internalModel.bounding_constr_mip[20][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[20][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[20][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[21][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[21][:vars] == Any[:(x[39])]
    @test m.internalModel.bounding_constr_mip[21][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[21][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[21][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[22][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[22][:vars] == Any[:(x[42])]
    @test m.internalModel.bounding_constr_mip[22][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[22][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[22][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[23][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[23][:vars] == Any[:(x[45])]
    @test m.internalModel.bounding_constr_mip[23][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[23][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[23][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[24][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[24][:vars] == Any[:(x[46])]
    @test m.internalModel.bounding_constr_mip[24][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[24][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[24][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[25][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[25][:vars] == Any[:(x[47])]
    @test m.internalModel.bounding_constr_mip[25][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[25][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[25][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[26][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[26][:vars] == Any[:(x[48])]
    @test m.internalModel.bounding_constr_mip[26][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[26][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[26][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[27][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[27][:vars] == Any[:(x[49])]
    @test m.internalModel.bounding_constr_mip[27][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[27][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[27][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[28][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[28][:vars] == Any[:(x[50])]
    @test m.internalModel.bounding_constr_mip[28][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[28][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[28][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[29][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[29][:vars] == Any[:(x[51])]
    @test m.internalModel.bounding_constr_mip[29][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[29][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[29][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[30][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[30][:vars] == Any[:(x[52])]
    @test m.internalModel.bounding_constr_mip[30][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[30][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[30][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[31][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[31][:vars] == Any[:(x[53])]
    @test m.internalModel.bounding_constr_mip[31][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[31][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[31][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[32][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[32][:vars] == Any[:(x[54])]
    @test m.internalModel.bounding_constr_mip[32][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[32][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[32][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[33][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[33][:vars] == Any[:(x[55])]
    @test m.internalModel.bounding_constr_mip[33][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[33][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[33][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[34][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[34][:vars] == Any[:(x[56])]
    @test m.internalModel.bounding_constr_mip[34][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[34][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[34][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[35][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[35][:vars] == Any[:(x[57])]
    @test m.internalModel.bounding_constr_mip[35][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[35][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[35][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[36][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[36][:vars] == Any[:(x[58])]
    @test m.internalModel.bounding_constr_mip[36][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[36][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[36][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[37][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[37][:vars] == Any[:(x[59])]
    @test m.internalModel.bounding_constr_mip[37][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[37][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[37][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[38][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[38][:vars] == Any[:(x[60])]
    @test m.internalModel.bounding_constr_mip[38][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[38][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[38][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[39][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[39][:vars] == Any[:(x[61])]
    @test m.internalModel.bounding_constr_mip[39][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[39][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[39][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[40][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[40][:vars] == Any[:(x[62])]
    @test m.internalModel.bounding_constr_mip[40][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[40][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[40][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[41][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[41][:vars] == Any[:(x[63])]
    @test m.internalModel.bounding_constr_mip[41][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[41][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[41][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[42][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[42][:vars] == Any[:(x[64])]
    @test m.internalModel.bounding_constr_mip[42][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[42][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[42][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[43][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[43][:vars] == Any[:(x[65])]
    @test m.internalModel.bounding_constr_mip[43][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[43][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[43][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[44][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[44][:vars] == Any[:(x[66])]
    @test m.internalModel.bounding_constr_mip[44][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[44][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[44][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[45][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[45][:vars] == Any[:(x[67])]
    @test m.internalModel.bounding_constr_mip[45][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[45][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[45][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[46][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[46][:vars] == Any[:(x[68])]
    @test m.internalModel.bounding_constr_mip[46][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[46][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[46][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[47][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[47][:vars] == Any[:(x[69])]
    @test m.internalModel.bounding_constr_mip[47][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[47][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[47][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[48][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[48][:vars] == Any[:(x[70])]
    @test m.internalModel.bounding_constr_mip[48][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[48][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[48][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[49][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[49][:vars] == Any[:(x[71])]
    @test m.internalModel.bounding_constr_mip[49][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[49][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[49][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[50][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[50][:vars] == Any[:(x[72])]
    @test m.internalModel.bounding_constr_mip[50][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[50][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[50][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[51][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[51][:vars] == Any[:(x[73])]
    @test m.internalModel.bounding_constr_mip[51][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[51][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[51][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[52][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[52][:vars] == Any[:(x[74])]
    @test m.internalModel.bounding_constr_mip[52][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[52][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[52][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[53][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[53][:vars] == Any[:(x[75])]
    @test m.internalModel.bounding_constr_mip[53][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[53][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[53][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[54][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[54][:vars] == Any[:(x[77])]
    @test m.internalModel.bounding_constr_mip[54][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[54][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[54][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[55][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[55][:vars] == Any[:(x[79])]
    @test m.internalModel.bounding_constr_mip[55][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[55][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[55][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[56][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[56][:vars] == Any[:(x[80])]
    @test m.internalModel.bounding_constr_mip[56][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[56][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[56][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[57][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[57][:vars] == Any[:(x[81])]
    @test m.internalModel.bounding_constr_mip[57][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[57][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[57][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[58][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[58][:vars] == Any[:(x[82])]
    @test m.internalModel.bounding_constr_mip[58][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[58][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[58][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[59][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[59][:vars] == Any[:(x[83])]
    @test m.internalModel.bounding_constr_mip[59][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[59][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[59][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[60][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[60][:vars] == Any[:(x[84])]
    @test m.internalModel.bounding_constr_mip[60][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[60][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[60][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[61][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[61][:vars] == Any[:(x[86])]
    @test m.internalModel.bounding_constr_mip[61][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[61][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[61][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[62][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[62][:vars] == Any[:(x[87])]
    @test m.internalModel.bounding_constr_mip[62][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[62][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[62][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[63][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[63][:vars] == Any[:(x[89])]
    @test m.internalModel.bounding_constr_mip[63][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[63][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[63][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[64][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[64][:vars] == Any[:(x[91])]
    @test m.internalModel.bounding_constr_mip[64][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[64][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[64][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[65][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[65][:vars] == Any[:(x[93])]
    @test m.internalModel.bounding_constr_mip[65][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[65][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[65][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[66][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[66][:vars] == Any[:(x[94])]
    @test m.internalModel.bounding_constr_mip[66][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[66][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[66][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[67][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[67][:vars] == Any[:(x[96])]
    @test m.internalModel.bounding_constr_mip[67][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[67][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[67][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[68][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[68][:vars] == Any[:(x[97])]
    @test m.internalModel.bounding_constr_mip[68][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[68][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[68][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[69][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[69][:vars] == Any[:(x[98])]
    @test m.internalModel.bounding_constr_mip[69][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[69][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[69][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[70][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[70][:vars] == Any[:(x[99])]
    @test m.internalModel.bounding_constr_mip[70][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[70][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[70][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[71][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[71][:vars] == Any[:(x[100])]
    @test m.internalModel.bounding_constr_mip[71][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[71][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[71][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[72][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[72][:vars] == Any[:(x[101])]
    @test m.internalModel.bounding_constr_mip[72][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[72][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[72][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[73][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[73][:vars] == Any[:(x[102])]
    @test m.internalModel.bounding_constr_mip[73][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[73][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[73][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[74][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[74][:vars] == Any[:(x[103])]
    @test m.internalModel.bounding_constr_mip[74][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[74][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[74][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[75][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[75][:vars] == Any[:(x[104])]
    @test m.internalModel.bounding_constr_mip[75][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[75][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[75][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[76][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[76][:vars] == Any[:(x[105])]
    @test m.internalModel.bounding_constr_mip[76][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[76][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[76][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[77][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[77][:vars] == Any[:(x[106])]
    @test m.internalModel.bounding_constr_mip[77][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[77][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[77][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[78][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[78][:vars] == Any[:(x[107])]
    @test m.internalModel.bounding_constr_mip[78][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[78][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[78][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[79][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[79][:vars] == Any[:(x[108])]
    @test m.internalModel.bounding_constr_mip[79][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[79][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[79][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[80][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[80][:vars] == Any[:(x[109])]
    @test m.internalModel.bounding_constr_mip[80][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[80][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[80][:cnt] == 1
end

@testset "Expression parsing corner cases" begin
    @testset "Corner Cases - 1 : sign convertor special case" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log=100)

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

    @testset "Corner Cases - 2 : full sub-expression" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(OutputFlag=0),
                                log=100)

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
