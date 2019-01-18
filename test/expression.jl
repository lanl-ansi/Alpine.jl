@testset "Expression Parsing || bilinear || Affine || exprs.jl" begin

    test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)

    m=exprstest(solver=test_solver)

    JuMP.build(m)

    ex = m.internalModel.bounding_constr_expr_mip[1]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-1.0]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[1][:coefs]
    @test affdict[:vars] == [:(x[1])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[1][:vars]
    @test isapprox(affdict[:rhs], 109.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[1][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[1][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[2]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [3.0,3.0,3.0,3.0]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[2][:coefs]
    @test affdict[:vars] == [:(x[8]),:(x[9]),:(x[10]),:(x[11])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[2][:vars]
    @test isapprox(affdict[:rhs], 111.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[2][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[2][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[3]
    affdict = Alpine.expr_linear_to_affine(ex)
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
    affdict = Alpine.expr_linear_to_affine(ex)
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
    affdict = Alpine.expr_linear_to_affine(ex)
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
    affdict = Alpine.expr_linear_to_affine(ex)
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
    affdict = Alpine.expr_linear_to_affine(ex)
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
    affdict = Alpine.expr_linear_to_affine(ex)
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
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [13.0,-1.0,30.0,1.0]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[9][:coefs]
    @test affdict[:vars] == [:(x[1]),:(x[2]),:(x[3]),:(x[4])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[9][:vars]
    @test isapprox(affdict[:rhs], 77.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[9][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[9][:sense]
end

@testset "Expression Parsing || bilinear || Affine || nlp1.jl" begin
    test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           loglevel=100)
    m=nlp1(solver=test_solver)

    JuMP.build(m)

    ex = m.internalModel.bounding_constr_expr_mip[1]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[1][:coefs]
    @test affdict[:vars] == [:(x[5])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[1][:vars]
    @test isapprox(affdict[:rhs], 8.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[1][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[1][:sense]
end

@testset "Expression Parsing || bilinear || Affine || nlp3.jl" begin

    test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)

    m=nlp3(solver=test_solver)

    JuMP.build(m)

    ex = m.internalModel.bounding_constr_expr_mip[1]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [0.0025,0.0025]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[1][:coefs]
    @test affdict[:vars] == [:(x[4]),:(x[6])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[1][:vars]
    @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[1][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[1][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[2]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [0.0025,-0.0025,0.0025]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[2][:coefs]
    @test affdict[:vars] == [:(x[5]),:(x[4]),:(x[7])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[2][:vars]
    @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[2][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[2][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[3]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [0.01, -0.01]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[3][:coefs]
    @test affdict[:vars] == [:(x[8]),:(x[5])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[3][:vars]
    @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[3][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[3][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[4]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test (affdict[:coefs] .== [100.0, -1.0, 833.33252]) == [true, true, true]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[4][:coefs]
    @test affdict[:vars] == [:(x[1]),:(x[9]),:(x[4])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[4][:vars]
    @test isapprox(affdict[:rhs], 83333.333; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[4][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[4][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[5]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0,-1.0,-1250.0,1250.0]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[5][:coefs]
    @test affdict[:vars] == [:(x[10]),:(x[11]),:(x[4]),:(x[5])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[5][:vars]
    @test isapprox(affdict[:rhs], 0.0; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[5][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[5][:sense]

    ex = m.internalModel.bounding_constr_expr_mip[6]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0,-1.0,-2500.0]
    @test affdict[:coefs] == m.internalModel.bounding_constr_mip[6][:coefs]
    @test affdict[:vars] == [:(x[12]),:(x[13]),:(x[5])]
    @test affdict[:vars] == m.internalModel.bounding_constr_mip[6][:vars]
    @test isapprox(affdict[:rhs], -1.25e6; atol = 1e-3)
    @test affdict[:rhs] == m.internalModel.bounding_constr_mip[6][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == m.internalModel.bounding_constr_mip[6][:sense]
end

@testset "Expression Parsing || bilinear || Simple || bi1.jl " begin

    test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           loglevel=100)

    m = operator_c(solver=test_solver)

    JuMP.build(m) # Setup internal model

    @test length(keys(m.internalModel.nonconvex_terms)) == 8
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 1), Expr(:ref, :x, 1)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 2)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 3)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 4), Expr(:ref, :x, 4)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 3)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 4)])

    # TODO setup detailed check on this problem
end

@testset "Expression Parsing || bilinear || Complex || blend029.jl " begin

    test_solver = AlpineSolver(minlp_solver=pavito_solver,mip_solver=CbcSolver(logLevel=0),loglevel=100)

    m = blend029(solver=test_solver)

    JuMP.build(m) # Setup internal model
    @test length(keys(m.internalModel.nonconvex_terms)) == 28
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 37), Expr(:ref, :x, 55)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 38), Expr(:ref, :x, 56)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 37), Expr(:ref, :x, 26)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 43), Expr(:ref, :x, 26)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 47), Expr(:ref, :x, 59)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 48), Expr(:ref, :x, 60)])
    @test haskey(m.internalModel.nonconvex_terms, [Expr(:ref, :x, 47), Expr(:ref, :x, 36)])

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

    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[55])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[55])]][:lifted_var_ref] == :(x[103])
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[55])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[55])]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_constr_mip[206][:rhs] == 0.0
    @test m.internalModel.bounding_constr_mip[206][:vars] == Any[:(x[107]), :(x[103]), :(x[108]), :(x[109]), :(x[110]), :(x[2]), :(x[14])]
    @test m.internalModel.bounding_constr_mip[206][:sense] == :(==)
    @test m.internalModel.bounding_constr_mip[206][:coefs] == Any[1.0, -1.0, 1.0, 1.0, 1.0, -0.6, -0.2]
    @test m.internalModel.bounding_constr_mip[206][:cnt] == 7

    @test m.internalModel.nonconvex_terms[[:(x[38]), :(x[56])]][:id] == 5
    @test m.internalModel.nonconvex_terms[[:(x[38]), :(x[56])]][:lifted_var_ref] == :(x[107])
    @test m.internalModel.nonconvex_terms[[:(x[38]), :(x[56])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[38]), :(x[56])]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[26])]][:id] == 6
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[26])]][:lifted_var_ref] == :(x[108])
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[26])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[26])]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[32])]][:id] == 8
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[32])]][:lifted_var_ref] == :(x[110])
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[32])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[37]), :(x[32])]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_constr_mip[213][:rhs] == 0.0
    @test m.internalModel.bounding_constr_mip[213][:vars] == Any[:(x[129]), :(x[127]), :(x[124]), :(x[130]), :(x[6]), :(x[18])]
    @test m.internalModel.bounding_constr_mip[213][:sense] == :(==)
    @test m.internalModel.bounding_constr_mip[213][:coefs] == Any[1.0, -1.0, -1.0, 1.0, -0.4, -0.4]
    @test m.internalModel.bounding_constr_mip[213][:cnt] == 6

    @test m.internalModel.nonconvex_terms[[:(x[48]), :(x[60])]][:id] == 27
    @test m.internalModel.nonconvex_terms[[:(x[48]), :(x[60])]][:lifted_var_ref] == :(x[129])
    @test m.internalModel.nonconvex_terms[[:(x[48]), :(x[60])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[48]), :(x[60])]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[59])]][:id] == 25
    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[59])]][:lifted_var_ref] == :(x[127])
    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[59])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[59])]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[36])]][:id] == 28
    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[36])]][:lifted_var_ref] == :(x[130])
    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[36])]][:convexified] == false
    @test m.internalModel.nonconvex_terms[[:(x[47]), :(x[36])]][:nonlinear_type] == :BILINEAR
end

@testset "Expression Parsing || multilinear || Simple || multi.jl " begin

    test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)

    m = multi3(solver=test_solver, exprmode=1)

    JuMP.build(m) # Setup internal model

    @test length(keys(m.internalModel.nonconvex_terms)) == 1
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:lifted_var_ref] == :(x[4])
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :MULTILINEAR

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

    @test length(keys(m.internalModel.nonconvex_terms)) == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[4])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[5])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi3(solver=test_solver, exprmode=3)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[4])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[5])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=1)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 1

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[5])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :MULTILINEAR

    @test m.internalModel.bounding_constr_mip[1][:rhs] == 4.0
    @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3]), :(x[4])]
    @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0]
    @test m.internalModel.bounding_constr_mip[1][:cnt] == 4

    m = multi4(solver=test_solver, exprmode=2)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 3

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:id] == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:nonlinear_type] == :BILINEAR


    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing


    m = multi4(solver=test_solver, exprmode=3)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 2

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[:(x[5]), :(x[3]), :(x[4])]][:id] == 2
    @test m.internalModel.nonconvex_terms[[:(x[5]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[:(x[5]), :(x[3]), :(x[4])]][:nonlinear_type] == :MULTILINEAR


    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=4)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 2

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[5])]][:id] == 2
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[5])]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[5])]][:nonlinear_type] == :MULTILINEAR

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=5)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 3

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[7])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=6)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 2

    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=7)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 3

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] == :BILINEAR


    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing


    m = multi4(solver=test_solver, exprmode=8)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 2

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[5]), :(x[4])]][:id] == 2
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[5]), :(x[4])]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[5]), :(x[4])]][:nonlinear_type] == :MULTILINEAR


    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing


    m = multi4(solver=test_solver, exprmode=9)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 2

    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] == :BILINEAR

    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=10)

    JuMP.build(m)
    @test length(keys(m.internalModel.nonconvex_terms)) == 3

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] == :(x[7])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] == :BILINEAR


    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing

    m = multi4(solver=test_solver, exprmode=11)

    JuMP.build(m)

    @test length(keys(m.internalModel.nonconvex_terms)) == 3

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] == :(x[6])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] == :(x[7])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] == :BILINEAR


    @test m.internalModel.bounding_obj_mip[:rhs] == 0
    @test m.internalModel.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test m.internalModel.bounding_obj_mip[:coefs] == [1.0]
    @test m.internalModel.bounding_obj_mip[:cnt] == 1
    @test m.internalModel.bounding_obj_mip[:sense] == nothing
end

@testset "Expression Parsing || bilinear || Complex-div || div.jl" begin
    test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)

    m = div(solver=test_solver)

    JuMP.build(m) # Setup internal model

    @test length(keys(m.internalModel.nonconvex_terms)) == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:id] == 1
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:lifted_var_ref] == :(x[3])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:nonlinear_type] == :MONOMIAL

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:id] == 2
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[4])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:nonlinear_type] == :MONOMIAL

    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 3
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] == :(x[5])
    @test m.internalModel.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] == :BILINEAR

    aff_mip = m.internalModel.bounding_constr_mip

    @test aff_mip[1][:rhs] == 0.0
    @test aff_mip[1][:vars] == Any[:(x[1])]
    @test aff_mip[1][:sense] == :(>=)
    @test round(aff_mip[1][:coefs][1]; digits=1) == -3.0
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

@testset "Expression Parsing || part1 " begin
    m = Model(solver=AlpineSolver(nlp_solver=IpoptSolver(),mip_solver=CbcSolver(logLevel=0),loglevel=100))
    @variable(m, x[1:4]>=0)
    @NLconstraint(m, x[1]^2 >= 1)  				  # Basic monomial x[5]=x[1]^2
    @NLconstraint(m, x[1]*x[2] <= 1)				  # x[6] <= 1 : x[6] = x[1]*x[2]
    @NLconstraint(m, x[1]^2 * x[2]^2 <= 1)          # x[5] + x[7] <= 1 : x[7] = x[2]^2
    @NLconstraint(m, x[1]*(x[2]*x[3]) >= 1)         # x[9] >= 1 : x[8] = x[2] * x[3] && x[9] = x[1]*x[8]
    @NLconstraint(m, x[1]^2*(x[2]^2 * x[3]^2) <= 1) # x[12] <= 1 : x[10] = x[3] ^ 2 && x[11] = x[7] * x[10] && x[12] = x[11]*[5]

    JuMP.build(m)

    @test m.internalModel.bounding_constr_expr_mip[1] == :(x[5]-1.0>=0.0)
    @test m.internalModel.bounding_constr_expr_mip[2] == :(x[6]-1.0<=0.0)
    @test m.internalModel.bounding_constr_expr_mip[3] == :(x[8]-1.0<=0.0)
    @test m.internalModel.bounding_constr_expr_mip[4] == :(x[10]-1.0>=0.0)
    @test m.internalModel.bounding_constr_expr_mip[5] == :(x[13]-1.0<=0.0)
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[1])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[2])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[2])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[3])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[5]), :(x[7])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[9])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[3])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[2])])
    @test haskey(m.internalModel.nonconvex_terms, [:(x[5]), :(x[12])])

    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2])]][:id] == 2
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 3
    @test m.internalModel.nonconvex_terms[[:(x[5]), :(x[7])]][:id] == 4
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[3])]][:id] == 5
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[9])]][:id] == 6
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 7
    @test m.internalModel.nonconvex_terms[[:(x[7]), :(x[11])]][:id] == 8
    @test m.internalModel.nonconvex_terms[[:(x[5]), :(x[12])]][:id] == 9
end

@testset "Expression Parsing || part2" begin
    m = Model(solver=AlpineSolver(nlp_solver=IpoptSolver(),mip_solver=CbcSolver(logLevel=0),loglevel=100))

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

    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[2])]) #5
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[5])]) #6
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[1])]) #7
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[2])]) #8
    @test haskey(m.internalModel.nonconvex_terms, [:(x[7]), :(x[8])]) #9
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[3])]) #10
    @test haskey(m.internalModel.nonconvex_terms, [:(x[9]), :(x[10])]) #11
    @test haskey(m.internalModel.nonconvex_terms, [:(x[8]), :(x[10])]) #12
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[12])]) #13
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[7])])  #14
    @test haskey(m.internalModel.nonconvex_terms, [:(x[14]), :(x[10])]) #15
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[3])]) #16
    @test haskey(m.internalModel.nonconvex_terms, [:(x[7]), :(x[16])]) #17
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[8])]) #18
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[18])]) #19

    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[5])]][:id] == 2
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 3
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 4
    @test m.internalModel.nonconvex_terms[[:(x[7]), :(x[8])]][:id] == 5
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 6
    @test m.internalModel.nonconvex_terms[[:(x[9]), :(x[10])]][:id] == 7
    @test m.internalModel.nonconvex_terms[[:(x[8]), :(x[10])]][:id] == 8
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[12])]][:id] == 9
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[7])]][:id] == 10
    @test m.internalModel.nonconvex_terms[[:(x[14]), :(x[10])]][:id] == 11
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[3])]][:id] == 12
    @test m.internalModel.nonconvex_terms[[:(x[7]), :(x[16])]][:id] == 13
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[8])]][:id] == 14
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[18])]][:id] == 15
end

@testset "Expression Parsing || part3" begin
    m = Model(solver=AlpineSolver(nlp_solver=IpoptSolver(),mip_solver=CbcSolver(logLevel=0),loglevel=100))

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

    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[2])]) #5
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[5])]) #6
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]), :(x[6])]) #7
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[1])]) #8
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[8])]) #9
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[9])]) #10
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]), :(x[10])]) #11
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[2])]) #12
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[12])])  #13
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[13])]) #14
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]), :(x[14])]) #15
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]), :(x[3])]) #16
    @test haskey(m.internalModel.nonconvex_terms, [:(x[5]), :(x[16])]) #17
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]), :(x[17])]) #18
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]), :(x[4])]) #19
    @test haskey(m.internalModel.nonconvex_terms, [:(x[6]), :(x[19])]) #20
    @test haskey(m.internalModel.nonconvex_terms, [:(x[8]), :(x[12])]) #21
    @test haskey(m.internalModel.nonconvex_terms, [:(x[21]), :(x[16])]) #22
    @test haskey(m.internalModel.nonconvex_terms, [:(x[22]), :(x[19])]) #23

    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[2])]][:id] == 1
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[5])]][:id] == 2
    @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[6])]][:id] == 3
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 4
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[8])]][:id] == 5
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[9])]][:id] == 6
    @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[10])]][:id] == 7
    @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 8
    @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[12])]][:id] == 9
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[13])]][:id] == 10
    @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[14])]][:id] == 11
    @test m.internalModel.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 12
    @test m.internalModel.nonconvex_terms[[:(x[5]), :(x[16])]][:id] == 13
    @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[17])]][:id] == 14
    @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[4])]][:id] == 15
    @test m.internalModel.nonconvex_terms[[:(x[6]), :(x[19])]][:id] == 16
    @test m.internalModel.nonconvex_terms[[:(x[8]), :(x[12])]][:id] == 17
    @test m.internalModel.nonconvex_terms[[:(x[21]), :(x[16])]][:id] == 18
    @test m.internalModel.nonconvex_terms[[:(x[22]), :(x[19])]][:id] == 19
end

@testset "Expression Parsing || part7" begin
    m = Model(solver=AlpineSolver(nlp_solver=IpoptSolver(),mip_solver=CbcSolver(logLevel=0),loglevel=100))
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

    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[2]),:(x[3]),:(x[4])]) #5
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[1])]) #6
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]),:(x[3]),:(x[4]),:(x[6])]) #7
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]),:(x[2])]) #8
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[3]),:(x[4]),:(x[8])]) #9
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]),:(x[3])]) #10
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]),:(x[4])]) #11
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[2]),:(x[10]),:(x[11])]) #12
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[4]),:(x[8]),:(x[10])]) #13
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]),:(x[3]),:(x[6]),:(x[11])]) #14
    @test haskey(m.internalModel.nonconvex_terms, [:(x[6]),:(x[8]),:(x[10]),:(x[11])]) #15

    @test m.internalModel.nonconvex_terms[[:(x[1]),:(x[2]),:(x[3]),:(x[4])]][:id] == 1 #5
    @test m.internalModel.nonconvex_terms[[:(x[1]),:(x[1])]][:id] == 2 #6
    @test m.internalModel.nonconvex_terms[[:(x[2]),:(x[3]),:(x[4]),:(x[6])]][:id] == 3 #7
    @test m.internalModel.nonconvex_terms[[:(x[2]),:(x[2])]][:id] == 4 #8
    @test m.internalModel.nonconvex_terms[[:(x[1]),:(x[3]),:(x[4]),:(x[8])]][:id] == 5 #9
    @test m.internalModel.nonconvex_terms[[:(x[3]),:(x[3])]][:id] == 6 #10
    @test m.internalModel.nonconvex_terms[[:(x[4]),:(x[4])]][:id] == 7 #11
    @test m.internalModel.nonconvex_terms[[:(x[1]),:(x[2]),:(x[10]),:(x[11])]][:id] ==  8  #12
    @test m.internalModel.nonconvex_terms[[:(x[1]),:(x[4]),:(x[8]),:(x[10])]][:id] == 9 #13
    @test m.internalModel.nonconvex_terms[[:(x[2]),:(x[3]),:(x[6]),:(x[11])]][:id] == 10 #14
    @test m.internalModel.nonconvex_terms[[:(x[6]),:(x[8]),:(x[10]),:(x[11])]][:id] == 11 #15
end

@testset "Expression Parsing || part8" begin
    m = Model(solver=AlpineSolver(nlp_solver=IpoptSolver(),
           mip_solver=CbcSolver(logLevel=0),
           loglevel=100))
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

    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[2]),:(x[3])]) #5
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]),:(x[5])]) #6
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[1])]) #7
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]),:(x[3]),:(x[7])]) #8
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]),:(x[8])]) #9
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]),:(x[2])]) #10
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]),:(x[10])]) #11
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[4]),:(x[11])]) #12
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]),:(x[3])]) #13
    @test haskey(m.internalModel.nonconvex_terms, [:(x[2]),:(x[13])]) #14
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[4]),:(x[14])]) #15
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[10])]) #16
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]),:(x[4]),:(x[16])]) #17
    @test haskey(m.internalModel.nonconvex_terms, [:(x[1]),:(x[2])]) #18
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]),:(x[18]),:(x[13])]) #19
    @test haskey(m.internalModel.nonconvex_terms, [:(x[4]),:(x[4])]) #20
    @test haskey(m.internalModel.nonconvex_terms, [:(x[3]),:(x[18]),:(x[20])]) #21
    @test haskey(m.internalModel.nonconvex_terms, [:(x[18]),:(x[13]),:(x[20])]) #22
    @test haskey(m.internalModel.nonconvex_terms, [:(x[7]),:(x[10]),:(x[13])]) #23
    @test haskey(m.internalModel.nonconvex_terms, [:(x[23]),:(x[20])]) #24
end

@testset "Expression Parsing || Convex" begin

    @testset "Convex Parsing :: PART I" begin

        test_solver = AlpineSolver(nlp_solver=IpoptSolver(),mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = convex_test(test_solver)

        JuMP.build(m)

        @test m.internalModel.num_constr_convex == 21

        # 0 : OBJ
        @test m.internalModel.obj_structure == :convex
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
        @test m.internalModel.constr_structure[1] == :convex
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
        @test m.internalModel.constr_structure[2] == :convex
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
        @test m.internalModel.constr_structure[4] == :convex
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
        @test m.internalModel.constr_structure[5] == :convex
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
        @test m.internalModel.constr_structure[6] == :convex
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
        @test m.internalModel.constr_structure[7] == :convex
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
        @test m.internalModel.constr_structure[8] == :convex
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
        @test m.internalModel.constr_structure[9] == :convex
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
        @test m.internalModel.constr_structure[11] == :convex
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
        @test m.internalModel.constr_structure[14] == :convex
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
        @test m.internalModel.constr_structure[15] == :convex
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
        @test m.internalModel.constr_structure[19] == :convex
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
        @test m.internalModel.constr_structure[22] == :convex
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
        @test m.internalModel.constr_structure[25] == :convex
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
        @test m.internalModel.constr_structure[26] == :convex
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

@testset "Expression Prasing || Linear Lifting" begin
    @testset "Expression Parsing || Linear Lifting || nlp2" begin
        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(),
                               disc_ratio=8,
                               loglevel=100)

        m = nlp2(solver=test_solver)

        JuMP.build(m)

        @test length(m.internalModel.linear_terms) == 2
        @test length(m.internalModel.nonconvex_terms) == 4

        lk = Vector{Any}(undef, 2)
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


        @test haskey(m.internalModel.nonconvex_terms, [:(x[2]), :(x[2])])
        @test haskey(m.internalModel.nonconvex_terms, [:(x[4]), :(x[4])])
        @test haskey(m.internalModel.nonconvex_terms, [:(x[7]), :(x[7])])
        @test haskey(m.internalModel.nonconvex_terms, [:(x[1]), :(x[1])])
        @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 1
        @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 3
        @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[4])]][:id] == 2
        @test m.internalModel.nonconvex_terms[[:(x[7]), :(x[7])]][:id] == 4
        @test m.internalModel.nonconvex_terms[[:(x[1]), :(x[1])]][:lifted_var_ref].args[2] == 3
        @test m.internalModel.nonconvex_terms[[:(x[2]), :(x[2])]][:lifted_var_ref].args[2] == 6
        @test m.internalModel.nonconvex_terms[[:(x[4]), :(x[4])]][:lifted_var_ref].args[2] == 5
        @test m.internalModel.nonconvex_terms[[:(x[7]), :(x[7])]][:lifted_var_ref].args[2] == 8
    end

    @testset "Expression Parsing || Linear Lifting || general" begin
        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(logLevel=0),
                               loglevel=100)

        m = basic_linear_lift(solver=test_solver)

        JuMP.build(m) # Setup internal model

        lk = Vector{Any}(undef, 5)
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

        @test m.internalModel.nonconvex_terms[nlk1][:id] == 7
        @test m.internalModel.nonconvex_terms[nlk2][:id] == 3
        @test m.internalModel.nonconvex_terms[nlk3][:id] == 4
        @test m.internalModel.nonconvex_terms[nlk4][:id] == 6
        @test m.internalModel.nonconvex_terms[nlk5][:id] == 2
        @test m.internalModel.nonconvex_terms[nlk6][:id] == 5
        @test m.internalModel.nonconvex_terms[nlk7][:id] == 1
        @test m.internalModel.nonconvex_terms[nlk8][:id] == 9
        @test m.internalModel.nonconvex_terms[nlk9][:id] == 8

        @test m.internalModel.nonconvex_terms[nlk1][:lifted_var_ref].args[2] == 13
        @test m.internalModel.nonconvex_terms[nlk2][:lifted_var_ref].args[2] == 6
        @test m.internalModel.nonconvex_terms[nlk3][:lifted_var_ref].args[2] == 7
        @test m.internalModel.nonconvex_terms[nlk4][:lifted_var_ref].args[2] == 11
        @test m.internalModel.nonconvex_terms[nlk5][:lifted_var_ref].args[2] == 5
        @test m.internalModel.nonconvex_terms[nlk6][:lifted_var_ref].args[2] == 10
        @test m.internalModel.nonconvex_terms[nlk7][:lifted_var_ref].args[2] == 4
        @test m.internalModel.nonconvex_terms[nlk8][:lifted_var_ref].args[2] == 17
        @test m.internalModel.nonconvex_terms[nlk9][:lifted_var_ref].args[2] == 15

        @test m.internalModel.nonconvex_terms[nlk1][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[nlk2][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[nlk3][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk4][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[nlk5][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk6][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk7][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk8][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk9][:nonlinear_type] == :MONOMIAL

        @test m.internalModel.nonconvex_terms[nlk1][:var_idxs] == [8,9,12]
        @test m.internalModel.nonconvex_terms[nlk2][:var_idxs] == [2,2]
        @test m.internalModel.nonconvex_terms[nlk3][:var_idxs] == [2,3]
        @test m.internalModel.nonconvex_terms[nlk4][:var_idxs] == [8]
        @test m.internalModel.nonconvex_terms[nlk5][:var_idxs] == [1,3]
        @test m.internalModel.nonconvex_terms[nlk6][:var_idxs] == [8,9]
        @test m.internalModel.nonconvex_terms[nlk7][:var_idxs] == [1,2]
        @test m.internalModel.nonconvex_terms[nlk8][:var_idxs] == [16, 15]
        @test m.internalModel.nonconvex_terms[nlk9][:var_idxs] == [14]
    end

    @testset "Expression Parsing || complex || Affine || operator_b" begin
        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                loglevel=100)

        m=operator_b(solver=test_solver)

        JuMP.build(m)

        nlk = Vector{Any}(undef, 31)
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

        @test length(keys(m.internalModel.nonconvex_terms)) == 31
        ids = [1:31;]
        yids = [8,9,11,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,37,38,39,41,43,45,47,49,50]
        nltypes = [:BILINEAR, :BILINEAR, :BILINEAR, :BILINEAR, :MONOMIAL,
                   :MONOMIAL, :MULTILINEAR, :MONOMIAL, :MULTILINEAR, :MULTILINEAR,
                   :MULTILINEAR, :BILINEAR, :BILINEAR, :MULTILINEAR, :BILINEAR,
                   :BILINEAR, :MULTILINEAR, :MULTILINEAR, :BILINEAR, :MONOMIAL,
                   :MULTILINEAR, :MULTILINEAR, :MULTILINEAR, :sin, :cos,
                   :sin, :cos, :sin, :sin, :cos, :sin]
        for i in 1:31
            @test m.internalModel.nonconvex_terms[nlk[i]][:id] == i
            @test m.internalModel.nonconvex_terms[nlk[i]][:y_idx] == yids[i]
            @test m.internalModel.nonconvex_terms[nlk[i]][:nonlinear_type] == nltypes[i]
            @test m.internalModel.nonconvex_terms[nlk[i]][:y_type] == :Cont
        end

        lk = Vector{Any}(undef, 12)
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

    @testset "Expression Parsing || Linear Lifting || brainpc3" begin

        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                            mip_solver=CbcSolver(logLevel=0),
                            disc_ratio=8,
                            loglevel=100)

        m = brainpc3(solver=test_solver)

        JuMP.build(m)

        @test m.internalModel.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:y_idx] == 6913
        @test m.internalModel.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:id] == 2
        @test m.internalModel.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:lifted_constr_ref] == :(x[6913] == (*)(x[6912]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:y_idx] == 6970
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:id] == 25
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:lifted_constr_ref] == :(x[6970] == x[6903] * x[6969])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:y_idx] == 6915
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:id] == 3
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:lifted_constr_ref] == :(x[6915] == x[6903] * x[6914])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:y_idx] == 6920
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:id] == 5
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:lifted_constr_ref] == :(x[6920] == x[6903] * x[6919])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:y_idx] == 6928
        @test m.internalModel.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:id] == 8
        @test m.internalModel.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:lifted_constr_ref] == :(x[6928] == (*)(x[6927]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:y_idx] == 6953
        @test m.internalModel.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:id] == 18
        @test m.internalModel.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:lifted_constr_ref] == :(x[6953] == (*)(x[6952]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:y_idx] == 6963
        @test m.internalModel.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:id] == 22
        @test m.internalModel.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:lifted_constr_ref] == :(x[6963] == (*)(x[6962]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:y_idx] == 6975
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:id] == 27
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:lifted_constr_ref] == :(x[6975] == x[6903] * x[6974])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:y_idx] == 6935
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:id] == 11
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:lifted_constr_ref] == :(x[6935] == x[6903] * x[6934])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:y_idx] == 6960
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:id] == 21
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:lifted_constr_ref] == :(x[6960] == x[6903] * x[6959])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:y_idx] == 7015
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:id] == 43
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:lifted_constr_ref] == :(x[7015] == x[6903] * x[7014])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:y_idx] == 6940
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:id] == 13
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:lifted_constr_ref] == :(x[6940] == x[6903] * x[6939])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:y_idx] == 7018
        @test m.internalModel.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:id] == 44
        @test m.internalModel.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:lifted_constr_ref] == :(x[7018] == (*)(x[7017]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:y_idx] == 6930
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:id] == 9
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:lifted_constr_ref] == :(x[6930] == x[6903] * x[6929])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:y_idx] == 7013
        @test m.internalModel.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:id] == 42
        @test m.internalModel.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:lifted_constr_ref] == :(x[7013] == (*)(x[7012]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:y_idx] == 6948
        @test m.internalModel.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:id] == 16
        @test m.internalModel.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:lifted_constr_ref] == :(x[6948] == (*)(x[6947]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:y_idx] == 7033
        @test m.internalModel.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:id] == 50
        @test m.internalModel.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:lifted_constr_ref] == :(x[7033] == (*)(x[7032]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:y_idx] == 6998
        @test m.internalModel.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:id] == 36
        @test m.internalModel.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:lifted_constr_ref] == :(x[6998] == (*)(x[6997]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:y_idx] == 7028
        @test m.internalModel.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:id] == 48
        @test m.internalModel.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:lifted_constr_ref] == :(x[7028] == (*)(x[7027]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:y_idx] == 7038
        @test m.internalModel.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:id] == 52
        @test m.internalModel.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:lifted_constr_ref] == :(x[7038] == (*)(x[7037]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:y_idx] == 7008
        @test m.internalModel.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:id] == 40
        @test m.internalModel.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:lifted_constr_ref] == :(x[7008] == (*)(x[7007]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:y_idx] == 6965
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:id] == 23
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:lifted_constr_ref] == :(x[6965] == x[6903] * x[6964])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:y_idx] == 6925
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:id] == 7
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:lifted_constr_ref] == :(x[6925] == x[6903] * x[6924])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:y_idx] == 6980
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:id] == 29
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:lifted_constr_ref] == :(x[6980] == x[6903] * x[6979])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:y_idx] == 6995
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:id] == 35
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:lifted_constr_ref] == :(x[6995] == x[6903] * x[6994])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:y_idx] == 6973
        @test m.internalModel.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:id] == 26
        @test m.internalModel.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:lifted_constr_ref] == :(x[6973] == (*)(x[6972]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:y_idx] == 7020
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:id] == 45
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:lifted_constr_ref] == :(x[7020] == x[6903] * x[7019])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:y_idx] == 6993
        @test m.internalModel.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:id] == 34
        @test m.internalModel.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:lifted_constr_ref] == :(x[6993] == (*)(x[6992]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:y_idx] == 7003
        @test m.internalModel.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:id] == 38
        @test m.internalModel.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:lifted_constr_ref] == :(x[7003] == (*)(x[7002]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:y_idx] == 7023
        @test m.internalModel.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:id] == 46
        @test m.internalModel.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:lifted_constr_ref] == :(x[7023] == (*)(x[7022]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:y_idx] == 6909
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:id] == 1
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:lifted_constr_ref] == :(x[6909] == x[6903] * x[6908])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:y_idx] == 6950
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:id] == 17
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:lifted_constr_ref] == :(x[6950] == x[6903] * x[6949])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:y_idx] == 7005
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:id] == 39
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:lifted_constr_ref] == :(x[7005] == x[6903] * x[7004])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:y_idx] == 6918
        @test m.internalModel.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:id] == 4
        @test m.internalModel.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:lifted_constr_ref] == :(x[6918] == (*)(x[6917]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:y_idx] == 6938
        @test m.internalModel.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:id] == 12
        @test m.internalModel.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:lifted_constr_ref] == :(x[6938] == (*)(x[6937]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:y_idx] == 6933
        @test m.internalModel.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:id] == 10
        @test m.internalModel.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:lifted_constr_ref] == :(x[6933] == (*)(x[6932]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:y_idx] == 6990
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:id] == 33
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:lifted_constr_ref] == :(x[6990] == x[6903] * x[6989])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:y_idx] == 6978
        @test m.internalModel.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:id] == 28
        @test m.internalModel.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:lifted_constr_ref] == :(x[6978] == (*)(x[6977]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:y_idx] == 6988
        @test m.internalModel.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:id] == 32
        @test m.internalModel.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:lifted_constr_ref] == :(x[6988] == (*)(x[6987]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:y_idx] == 7035
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:id] == 51
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:lifted_constr_ref] == :(x[7035] == x[6903] * x[7034])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:y_idx] == 6968
        @test m.internalModel.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:id] == 24
        @test m.internalModel.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:lifted_constr_ref] == :(x[6968] == (*)(x[6967]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:y_idx] == 6958
        @test m.internalModel.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:id] == 20
        @test m.internalModel.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:lifted_constr_ref] == :(x[6958] == (*)(x[6957]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:y_idx] == 7010
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:id] == 41
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:lifted_constr_ref] == :(x[7010] == x[6903] * x[7009])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:y_idx] == 6943
        @test m.internalModel.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:id] == 14
        @test m.internalModel.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:lifted_constr_ref] == :(x[6943] == (*)(x[6942]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:y_idx] == 6983
        @test m.internalModel.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:id] == 30
        @test m.internalModel.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:lifted_constr_ref] == :(x[6983] == (*)(x[6982]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:y_idx] == 7000
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:id] == 37
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:lifted_constr_ref] == :(x[7000] == x[6903] * x[6999])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:y_idx] == 6945
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:id] == 15
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:lifted_constr_ref] == :(x[6945] == x[6903] * x[6944])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:y_idx] == 7025
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:id] == 47
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:lifted_constr_ref] == :(x[7025] == x[6903] * x[7024])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:y_idx] == 6955
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:id] == 19
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:lifted_constr_ref] == :(x[6955] == x[6903] * x[6954])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:y_idx] == 6923
        @test m.internalModel.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:id] == 6
        @test m.internalModel.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:lifted_constr_ref] == :(x[6923] == (*)(x[6922]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:y_idx] == 6985
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:id] == 31
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:lifted_constr_ref] == :(x[6985] == x[6903] * x[6984])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:y_idx] == 7030
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:id] == 49
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:lifted_constr_ref] == :(x[7030] == x[6903] * x[7029])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:constr_id] == Set(Any[0])
        lk1 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6910), (-1.0, 6909)])))
        @test m.internalModel.linear_terms[lk1][:y_idx] == 6911
        @test m.internalModel.linear_terms[lk1][:id] == 3
        @test m.internalModel.linear_terms[lk1][:y_type] == :(Cont)
        lk2 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3467), (1.0, 16)])))
        @test m.internalModel.linear_terms[lk2][:y_idx] == 6908
        @test m.internalModel.linear_terms[lk2][:id] == 1
        @test m.internalModel.linear_terms[lk2][:y_type] == :(Cont)
        lk3 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.001548, 6903), (1.0, 7016), (1.0, 5702)])))
        @test m.internalModel.linear_terms[lk3][:y_idx] == 7017
        @test m.internalModel.linear_terms[lk3][:id] == 67
        @test m.internalModel.linear_terms[lk3][:y_type] == :(Cont)
        lk4 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 5162), (1.0, 1711)])))
        @test m.internalModel.linear_terms[lk4][:y_idx] == 7004
        @test m.internalModel.linear_terms[lk4][:id] == 59
        @test m.internalModel.linear_terms[lk4][:y_type] == :(Cont)
        lk5 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7021), (0.001435, 6903), (1.0, 6002)])))
        @test m.internalModel.linear_terms[lk5][:y_idx] == 7022
        @test m.internalModel.linear_terms[lk5][:id] == 70
        @test m.internalModel.linear_terms[lk5][:y_type] == :(Cont)
        lk6 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6935), (1.0, 166)])))
        @test m.internalModel.linear_terms[lk6][:y_idx] == 6936
        @test m.internalModel.linear_terms[lk6][:id] == 18
        @test m.internalModel.linear_terms[lk6][:y_type] == :(Cont)
        lk7 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4262), (0.00247, 6903), (1.0, 6981)])))
        @test m.internalModel.linear_terms[lk7][:y_idx] == 6982
        @test m.internalModel.linear_terms[lk7][:id] == 46
        @test m.internalModel.linear_terms[lk7][:y_type] == :(Cont)
        lk8 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7005), (1.0, 1711)])))
        @test m.internalModel.linear_terms[lk8][:y_idx] == 7006
        @test m.internalModel.linear_terms[lk8][:id] == 60
        @test m.internalModel.linear_terms[lk8][:y_type] == :(Cont)
        lk9 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1951), (1.0, 5402)])))
        @test m.internalModel.linear_terms[lk9][:y_idx] == 7009
        @test m.internalModel.linear_terms[lk9][:id] == 62
        @test m.internalModel.linear_terms[lk9][:y_type] == :(Cont)
        lk10 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 226), (-1.0, 6945)])))
        @test m.internalModel.linear_terms[lk10][:y_idx] == 6946
        @test m.internalModel.linear_terms[lk10][:id] == 24
        @test m.internalModel.linear_terms[lk10][:y_type] == :(Cont)
        lk11 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7030), (1.0, 3151)])))
        @test m.internalModel.linear_terms[lk11][:y_idx] == 7031
        @test m.internalModel.linear_terms[lk11][:id] == 75
        @test m.internalModel.linear_terms[lk11][:y_type] == :(Cont)
        lk12 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6991), (1.0, 4622), (0.002066, 6903)])))
        @test m.internalModel.linear_terms[lk12][:y_idx] == 6992
        @test m.internalModel.linear_terms[lk12][:id] == 52
        @test m.internalModel.linear_terms[lk12][:y_type] == :(Cont)
        lk13 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7026), (1.0, 6302), (0.001323, 6903)])))
        @test m.internalModel.linear_terms[lk13][:y_idx] == 7027
        @test m.internalModel.linear_terms[lk13][:id] == 73
        @test m.internalModel.linear_terms[lk13][:y_type] == :(Cont)
        lk14 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 46), (1.0, 3497)])))
        @test m.internalModel.linear_terms[lk14][:y_idx] == 6914
        @test m.internalModel.linear_terms[lk14][:id] == 5
        @test m.internalModel.linear_terms[lk14][:y_type] == :(Cont)
        lk15 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6960), (1.0, 391)])))
        @test m.internalModel.linear_terms[lk15][:y_idx] == 6961
        @test m.internalModel.linear_terms[lk15][:id] == 33
        @test m.internalModel.linear_terms[lk15][:y_type] == :(Cont)
        lk16 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 631), (-1.0, 6975)])))
        @test m.internalModel.linear_terms[lk16][:y_idx] == 6976
        @test m.internalModel.linear_terms[lk16][:id] == 42
        @test m.internalModel.linear_terms[lk16][:y_type] == :(Cont)
        lk17 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1351), (-1.0, 6995)])))
        @test m.internalModel.linear_terms[lk17][:y_idx] == 6996
        @test m.internalModel.linear_terms[lk17][:id] == 54
        @test m.internalModel.linear_terms[lk17][:y_type] == :(Cont)
        lk18 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3557), (1.0, 106)])))
        @test m.internalModel.linear_terms[lk18][:y_idx] == 6924
        @test m.internalModel.linear_terms[lk18][:id] == 11
        @test m.internalModel.linear_terms[lk18][:y_type] == :(Cont)
        lk19 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3647), (1.0, 6941), (0.004795, 6903)])))
        @test m.internalModel.linear_terms[lk19][:y_idx] == 6942
        @test m.internalModel.linear_terms[lk19][:id] == 22
        @test m.internalModel.linear_terms[lk19][:y_type] == :(Cont)
        lk20 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3722), (1.0, 271)])))
        @test m.internalModel.linear_terms[lk20][:y_idx] == 6949
        @test m.internalModel.linear_terms[lk20][:id] == 26
        @test m.internalModel.linear_terms[lk20][:y_type] == :(Cont)
        lk21 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3647), (1.0, 196)])))
        @test m.internalModel.linear_terms[lk21][:y_idx] == 6939
        @test m.internalModel.linear_terms[lk21][:id] == 20
        @test m.internalModel.linear_terms[lk21][:y_type] == :(Cont)
        lk22 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3467), (1.0, 6911), (6.34e-6, 6903)])))
        @test m.internalModel.linear_terms[lk22][:y_idx] == 6912
        @test m.internalModel.linear_terms[lk22][:id] == 4
        @test m.internalModel.linear_terms[lk22][:y_type] == :(Cont)
        lk23 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1531), (1.0, 4982)])))
        @test m.internalModel.linear_terms[lk23][:y_idx] == 6999
        @test m.internalModel.linear_terms[lk23][:id] == 56
        @test m.internalModel.linear_terms[lk23][:y_type] == :(Cont)
        lk24 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 106), (-1.0, 6925)])))
        @test m.internalModel.linear_terms[lk24][:y_idx] == 6926
        @test m.internalModel.linear_terms[lk24][:id] == 12
        @test m.internalModel.linear_terms[lk24][:y_type] == :(Cont)
        lk25 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2251), (-1.0, 7015)])))
        @test m.internalModel.linear_terms[lk25][:y_idx] == 7016
        @test m.internalModel.linear_terms[lk25][:id] == 66
        @test m.internalModel.linear_terms[lk25][:y_type] == :(Cont)
        lk26 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6936), (0.004985, 6903), (1.0, 3617)])))
        @test m.internalModel.linear_terms[lk26][:y_idx] == 6937
        @test m.internalModel.linear_terms[lk26][:id] == 19
        @test m.internalModel.linear_terms[lk26][:y_type] == :(Cont)
        lk27 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6915), (1.0, 46)])))
        @test m.internalModel.linear_terms[lk27][:y_idx] == 6916
        @test m.internalModel.linear_terms[lk27][:id] == 6
        @test m.internalModel.linear_terms[lk27][:y_type] == :(Cont)
        lk28 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3557), (0.005968, 6903), (1.0, 6926)])))
        @test m.internalModel.linear_terms[lk28][:y_idx] == 6927
        @test m.internalModel.linear_terms[lk28][:id] == 13
        @test m.internalModel.linear_terms[lk28][:y_type] == :(Cont)
        lk29 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 271), (-1.0, 6950)])))
        @test m.internalModel.linear_terms[lk29][:y_idx] == 6951
        @test m.internalModel.linear_terms[lk29][:id] == 27
        @test m.internalModel.linear_terms[lk29][:y_type] == :(Cont)
        lk30 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.002971, 6903), (1.0, 6971), (1.0, 3962)])))
        @test m.internalModel.linear_terms[lk30][:y_idx] == 6972
        @test m.internalModel.linear_terms[lk30][:id] == 40
        @test m.internalModel.linear_terms[lk30][:y_type] == :(Cont)
        lk31 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 511), (-1.0, 6970)])))
        @test m.internalModel.linear_terms[lk31][:y_idx] == 6971
        @test m.internalModel.linear_terms[lk31][:id] == 39
        @test m.internalModel.linear_terms[lk31][:y_type] == :(Cont)
        lk32 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7031), (1.0, 6602), (0.00121, 6903)])))
        @test m.internalModel.linear_terms[lk32][:y_idx] == 7032
        @test m.internalModel.linear_terms[lk32][:id] == 76
        @test m.internalModel.linear_terms[lk32][:y_type] == :(Cont)
        lk33 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.00166, 6903), (1.0, 7011), (1.0, 5402)])))
        @test m.internalModel.linear_terms[lk33][:y_idx] == 7012
        @test m.internalModel.linear_terms[lk33][:id] == 64
        @test m.internalModel.linear_terms[lk33][:y_type] == :(Cont)
        lk34 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, -0.003214),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 16)])))
        @test m.internalModel.linear_terms[lk34][:y_idx] == 6910
        @test m.internalModel.linear_terms[lk34][:id] == 2
        @test m.internalModel.linear_terms[lk34][:y_type] == :(Cont)
        lk35 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7035), (1.0, 3451)])))
        @test m.internalModel.linear_terms[lk35][:y_idx] == 7036
        @test m.internalModel.linear_terms[lk35][:id] == 78
        @test m.internalModel.linear_terms[lk35][:y_type] == :(Cont)
        lk36 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2551), (1.0, 6002)])))
        @test m.internalModel.linear_terms[lk36][:y_idx] == 7019
        @test m.internalModel.linear_terms[lk36][:id] == 68
        @test m.internalModel.linear_terms[lk36][:y_type] == :(Cont)
        lk37 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3842), (1.0, 391)])))
        @test m.internalModel.linear_terms[lk37][:y_idx] == 6959
        @test m.internalModel.linear_terms[lk37][:id] == 32
        @test m.internalModel.linear_terms[lk37][:y_type] == :(Cont)
        lk38 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6966), (0.00307, 6903), (1.0, 3902)])))
        @test m.internalModel.linear_terms[lk38][:y_idx] == 6967
        @test m.internalModel.linear_terms[lk38][:id] == 37
        @test m.internalModel.linear_terms[lk38][:y_type] == :(Cont)
        lk39 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3722), (0.003829, 6903), (1.0, 6951)])))
        @test m.internalModel.linear_terms[lk39][:y_idx] == 6952
        @test m.internalModel.linear_terms[lk39][:id] == 28
        @test m.internalModel.linear_terms[lk39][:y_type] == :(Cont)
        lk40 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 451), (1.0, 3902)])))
        @test m.internalModel.linear_terms[lk40][:y_idx] == 6964
        @test m.internalModel.linear_terms[lk40][:id] == 35
        @test m.internalModel.linear_terms[lk40][:y_type] == :(Cont)
        lk41 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7006), (1.0, 5162), (0.001713, 6903)])))
        @test m.internalModel.linear_terms[lk41][:y_idx] == 7007
        @test m.internalModel.linear_terms[lk41][:id] == 61
        @test m.internalModel.linear_terms[lk41][:y_type] == :(Cont)
        lk42 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 811), (-1.0, 6980)])))
        @test m.internalModel.linear_terms[lk42][:y_idx] == 6981
        @test m.internalModel.linear_terms[lk42][:id] == 45
        @test m.internalModel.linear_terms[lk42][:y_type] == :(Cont)
        lk43 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7001), (1.0, 4982), (0.00173, 6903)])))
        @test m.internalModel.linear_terms[lk43][:y_idx] == 7002
        @test m.internalModel.linear_terms[lk43][:id] == 58
        @test m.internalModel.linear_terms[lk43][:y_type] == :(Cont)
        lk44 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6920), (1.0, 76)])))
        @test m.internalModel.linear_terms[lk44][:y_idx] == 6921
        @test m.internalModel.linear_terms[lk44][:id] == 9
        @test m.internalModel.linear_terms[lk44][:y_type] == :(Cont)
        lk45 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4622), (1.0, 1171)])))
        @test m.internalModel.linear_terms[lk45][:y_idx] == 6989
        @test m.internalModel.linear_terms[lk45][:id] == 50
        @test m.internalModel.linear_terms[lk45][:y_type] == :(Cont)
        lk46 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 991), (1.0, 4442)])))
        @test m.internalModel.linear_terms[lk46][:y_idx] == 6984
        @test m.internalModel.linear_terms[lk46][:id] == 47
        @test m.internalModel.linear_terms[lk46][:y_type] == :(Cont)
        lk47 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6930), (1.0, 136)])))
        @test m.internalModel.linear_terms[lk47][:y_idx] == 6931
        @test m.internalModel.linear_terms[lk47][:id] == 15
        @test m.internalModel.linear_terms[lk47][:y_type] == :(Cont)
        lk48 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4802), (0.001891, 6903), (1.0, 6996)])))
        @test m.internalModel.linear_terms[lk48][:y_idx] == 6997
        @test m.internalModel.linear_terms[lk48][:id] == 55
        @test m.internalModel.linear_terms[lk48][:y_type] == :(Cont)
        lk49 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 511), (1.0, 3962)])))
        @test m.internalModel.linear_terms[lk49][:y_idx] == 6969
        @test m.internalModel.linear_terms[lk49][:id] == 38
        @test m.internalModel.linear_terms[lk49][:y_type] == :(Cont)
        lk50 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 811), (1.0, 4262)])))
        @test m.internalModel.linear_terms[lk50][:y_idx] == 6979
        @test m.internalModel.linear_terms[lk50][:id] == 44
        @test m.internalModel.linear_terms[lk50][:y_type] == :(Cont)
        lk51 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.01551, 6903), (1.0, 6916), (1.0, 3497)])))
        @test m.internalModel.linear_terms[lk51][:y_idx] == 6917
        @test m.internalModel.linear_terms[lk51][:id] == 7
        @test m.internalModel.linear_terms[lk51][:y_type] == :(Cont)
        lk52 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 451), (-1.0, 6965)])))
        @test m.internalModel.linear_terms[lk52][:y_idx] == 6966
        @test m.internalModel.linear_terms[lk52][:id] == 36
        @test m.internalModel.linear_terms[lk52][:y_type] == :(Cont)
        lk53 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.005773, 6903), (1.0, 3587), (1.0, 6931)])))
        @test m.internalModel.linear_terms[lk53][:y_idx] == 6932
        @test m.internalModel.linear_terms[lk53][:id] == 16
        @test m.internalModel.linear_terms[lk53][:y_type] == :(Cont)
        lk54 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 631), (1.0, 4082)])))
        @test m.internalModel.linear_terms[lk54][:y_idx] == 6974
        @test m.internalModel.linear_terms[lk54][:id] == 41
        @test m.internalModel.linear_terms[lk54][:y_type] == :(Cont)
        lk55 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 991), (-1.0, 6985)])))
        @test m.internalModel.linear_terms[lk55][:y_idx] == 6986
        @test m.internalModel.linear_terms[lk55][:id] == 48
        @test m.internalModel.linear_terms[lk55][:y_type] == :(Cont)
        lk56 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4802), (1.0, 1351)])))
        @test m.internalModel.linear_terms[lk56][:y_idx] == 6994
        @test m.internalModel.linear_terms[lk56][:id] == 53
        @test m.internalModel.linear_terms[lk56][:y_type] == :(Cont)
        lk57 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6946), (0.004409, 6903), (1.0, 3677)])))
        @test m.internalModel.linear_terms[lk57][:y_idx] == 6947
        @test m.internalModel.linear_terms[lk57][:id] == 25
        @test m.internalModel.linear_terms[lk57][:y_type] == :(Cont)
        lk58 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2251), (1.0, 5702)])))
        @test m.internalModel.linear_terms[lk58][:y_idx] == 7014
        @test m.internalModel.linear_terms[lk58][:id] == 65
        @test m.internalModel.linear_terms[lk58][:y_type] == :(Cont)
        lk59 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7020), (1.0, 2551)])))
        @test m.internalModel.linear_terms[lk59][:y_idx] == 7021
        @test m.internalModel.linear_terms[lk59][:id] == 69
        @test m.internalModel.linear_terms[lk59][:y_type] == :(Cont)
        lk60 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6940), (1.0, 196)])))
        @test m.internalModel.linear_terms[lk60][:y_idx] == 6941
        @test m.internalModel.linear_terms[lk60][:id] == 21
        @test m.internalModel.linear_terms[lk60][:y_type] == :(Cont)
        lk61 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3527), (1.0, 6921), (0.008919, 6903)])))
        @test m.internalModel.linear_terms[lk61][:y_idx] == 6922
        @test m.internalModel.linear_terms[lk61][:id] == 10
        @test m.internalModel.linear_terms[lk61][:y_type] == :(Cont)
        lk62 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7025), (1.0, 2851)])))
        @test m.internalModel.linear_terms[lk62][:y_idx] == 7026
        @test m.internalModel.linear_terms[lk62][:id] == 72
        @test m.internalModel.linear_terms[lk62][:y_type] == :(Cont)
        lk63 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6986), (1.0, 4442), (0.002252, 6903)])))
        @test m.internalModel.linear_terms[lk63][:y_idx] == 6987
        @test m.internalModel.linear_terms[lk63][:id] == 49
        @test m.internalModel.linear_terms[lk63][:y_type] == :(Cont)
        lk64 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(0.002765, 6903), (1.0, 6976), (1.0, 4082)])))
        @test m.internalModel.linear_terms[lk64][:y_idx] == 6977
        @test m.internalModel.linear_terms[lk64][:id] == 43
        @test m.internalModel.linear_terms[lk64][:y_type] == :(Cont)
        lk65 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6956), (0.003269, 6903), (1.0, 3782)])))
        @test m.internalModel.linear_terms[lk65][:y_idx] == 6957
        @test m.internalModel.linear_terms[lk65][:id] == 31
        @test m.internalModel.linear_terms[lk65][:y_type] == :(Cont)
        lk66 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3587), (1.0, 136)])))
        @test m.internalModel.linear_terms[lk66][:y_idx] == 6929
        @test m.internalModel.linear_terms[lk66][:id] == 14
        @test m.internalModel.linear_terms[lk66][:y_type] == :(Cont)
        lk67 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 166), (1.0, 3617)])))
        @test m.internalModel.linear_terms[lk67][:y_idx] == 6934
        @test m.internalModel.linear_terms[lk67][:id] == 17
        @test m.internalModel.linear_terms[lk67][:y_type] == :(Cont)
        lk68 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7010), (1.0, 1951)])))
        @test m.internalModel.linear_terms[lk68][:y_idx] == 7011
        @test m.internalModel.linear_terms[lk68][:id] == 63
        @test m.internalModel.linear_terms[lk68][:y_type] == :(Cont)
        lk69 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7036), (0.001098, 6903), (1.0, 6902)])))
        @test m.internalModel.linear_terms[lk69][:y_idx] == 7037
        @test m.internalModel.linear_terms[lk69][:id] == 79
        @test m.internalModel.linear_terms[lk69][:y_type] == :(Cont)
        lk70 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3527), (1.0, 76)])))
        @test m.internalModel.linear_terms[lk70][:y_idx] == 6919
        @test m.internalModel.linear_terms[lk70][:id] == 8
        @test m.internalModel.linear_terms[lk70][:y_type] == :(Cont)
        lk71 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6961), (1.0, 3842), (0.00317, 6903)])))
        @test m.internalModel.linear_terms[lk71][:y_idx] == 6962
        @test m.internalModel.linear_terms[lk71][:id] == 34
        @test m.internalModel.linear_terms[lk71][:y_type] == :(Cont)
        lk72 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6990), (1.0, 1171)])))
        @test m.internalModel.linear_terms[lk72][:y_idx] == 6991
        @test m.internalModel.linear_terms[lk72][:id] == 51
        @test m.internalModel.linear_terms[lk72][:y_type] == :(Cont)
        lk73 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3451), (1.0, 6902)])))
        @test m.internalModel.linear_terms[lk73][:y_idx] == 7034
        @test m.internalModel.linear_terms[lk73][:id] == 77
        @test m.internalModel.linear_terms[lk73][:y_type] == :(Cont)
        lk74 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6302), (1.0, 2851)])))
        @test m.internalModel.linear_terms[lk74][:y_idx] == 7024
        @test m.internalModel.linear_terms[lk74][:id] == 71
        @test m.internalModel.linear_terms[lk74][:y_type] == :(Cont)
        lk75 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1531), (-1.0, 7000)])))
        @test m.internalModel.linear_terms[lk75][:y_idx] == 7001
        @test m.internalModel.linear_terms[lk75][:id] == 57
        @test m.internalModel.linear_terms[lk75][:y_type] == :(Cont)
        lk76 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3151), (1.0, 6602)])))
        @test m.internalModel.linear_terms[lk76][:y_idx] == 7029
        @test m.internalModel.linear_terms[lk76][:id] == 74
        @test m.internalModel.linear_terms[lk76][:y_type] == :(Cont)
        lk77 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3782), (1.0, 331)])))
        @test m.internalModel.linear_terms[lk77][:y_idx] == 6954
        @test m.internalModel.linear_terms[lk77][:id] == 29
        @test m.internalModel.linear_terms[lk77][:y_type] == :(Cont)
        lk78 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6955), (1.0, 331)])))
        @test m.internalModel.linear_terms[lk78][:y_idx] == 6956
        @test m.internalModel.linear_terms[lk78][:id] == 30
        @test m.internalModel.linear_terms[lk78][:y_type] == :(Cont)
        lk79 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 226), (1.0, 3677)])))
        @test m.internalModel.linear_terms[lk79][:y_idx] == 6944
        @test m.internalModel.linear_terms[lk79][:id] == 23
        @test m.internalModel.linear_terms[lk79][:y_type] == :(Cont)
    end
end

@testset "Expression Parsing || Basic Multiplication Operators (Machine Generated for diffs)" begin

    test_solver=AlpineSolver(nlp_solver=IpoptSolver(),
                           mip_solver=CbcSolver(logLevel=0),
                           loglevel=100)

    m = operator_basic(solver=test_solver)
    JuMP.build(m)

    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:y_idx] == 31
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:id] == 27
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[31] == x[3] * x[4])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:y_idx] == 83
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:id] == 79
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:lifted_constr_ref] == :(x[83] == x[4] * x[5] * x[76])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:y_idx] == 9
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:id] == 5
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:lifted_constr_ref] == :(x[9] == (*)(x[2]))
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:y_idx] == 19
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:id] == 15
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:lifted_constr_ref] == :(x[19] == x[5] * x[7])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:y_idx] == 32
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:id] == 28
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:lifted_constr_ref] == :(x[32] == x[2] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:y_idx] == 49
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:id] == 45
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:lifted_constr_ref] == :(x[49] == x[1] * x[2] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:y_idx] == 50
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:id] == 46
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:lifted_constr_ref] == :(x[50] == x[3] * x[5] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:y_idx] == 69
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:id] == 65
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:lifted_constr_ref] == :(x[69] == x[4] * x[47])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:y_idx] == 28
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:id] == 24
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:lifted_constr_ref] == :(x[28] == (*)(x[4]))
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:y_idx] == 37
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:id] == 33
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:lifted_constr_ref] == :(x[37] == x[4] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:y_idx] == 5
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:id] == 1
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:lifted_constr_ref] == :(x[5] == (*)(x[1]))
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:y_idx] == 96
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:id] == 92
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:lifted_constr_ref] == :(x[96] == x[5] * x[95])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:y_idx] == 41
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:id] == 37
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:lifted_constr_ref] == :(x[41] == x[2] * x[40])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:y_idx] == 26
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:id] == 22
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:lifted_constr_ref] == :(x[26] == x[6] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:y_idx] == 84
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:id] == 80
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:lifted_constr_ref] == :(x[84] == x[1] * x[78] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:y_idx] == 58
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:id] == 54
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:lifted_constr_ref] == :(x[58] == x[6] * x[43])
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:y_idx] == 106
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:id] == 102
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:lifted_constr_ref] == :(x[106] == x[4] * x[6] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:y_idx] == 90
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:id] == 86
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:lifted_constr_ref] == :(x[90] == x[2] * x[4] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:y_idx] == 42
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:id] == 38
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:lifted_constr_ref] == :(x[42] == x[1] * x[41])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:y_idx] == 38
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:id] == 34
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:lifted_constr_ref] == :(x[38] == x[2] * x[37])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:y_idx] == 87
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:id] == 83
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:lifted_constr_ref] == :(x[87] == x[5] * x[85])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:y_idx] == 66
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:id] == 62
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:lifted_constr_ref] == :(x[66] == x[2] * x[3] * x[5] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:y_idx] == 53
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:id] == 49
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[53] == x[5] * x[9] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:y_idx] == 107
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:id] == 103
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:lifted_constr_ref] == :(x[107] == x[3] * x[6] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:y_idx] == 17
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:id] == 13
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:lifted_constr_ref] == :(x[17] == x[2] * x[5])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:y_idx] == 57
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:id] == 53
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:lifted_constr_ref] == :(x[57] == x[6] * x[37])
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:y_idx] == 81
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:id] == 77
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:lifted_constr_ref] == :(x[81] == x[5] * x[7] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:y_idx] == 39
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:id] == 35
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:lifted_constr_ref] == :(x[39] == x[1] * x[38])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:y_idx] == 48
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:id] == 44
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:lifted_constr_ref] == :(x[48] == x[1] * x[3] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:y_idx] == 24
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:id] == 20
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:lifted_constr_ref] == :(x[24] == x[4] * x[23])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:y_idx] == 23
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:id] == 19
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:lifted_constr_ref] == :(x[23] == x[3] * x[17])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:y_idx] == 51
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:id] == 47
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[51] == x[1] * x[9] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:y_idx] == 75
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:id] == 71
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:lifted_constr_ref] == :(x[75] == x[4] * x[5] * x[7])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:y_idx] == 95
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:id] == 91
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[95] == x[9] * x[10] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:y_idx] == 11
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:id] == 7
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[11] == x[9] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:y_idx] == 100
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:id] == 96
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:lifted_constr_ref] == :(x[100] == x[1] * x[2] * x[37])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:y_idx] == 62
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:id] == 58
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:lifted_constr_ref] == :(x[62] == x[2] * x[3] * x[4] * x[5])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:y_idx] == 65
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:id] == 61
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[65] == x[1] * x[4] * x[9] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:y_idx] == 80
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:id] == 76
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:lifted_constr_ref] == :(x[80] == x[1] * x[7] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:y_idx] == 46
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:id] == 42
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:lifted_constr_ref] == :(x[46] == x[1] * x[2] * x[3])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:y_idx] == 59
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:id] == 55
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:lifted_constr_ref] == :(x[59] == x[17] * x[40])
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:y_idx] == 86
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:id] == 82
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:lifted_constr_ref] == :(x[86] == x[1] * x[85])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:y_idx] == 97
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:id] == 93
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:lifted_constr_ref] == :(x[97] == x[1] * x[2] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:y_idx] == 29
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:id] == 25
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:lifted_constr_ref] == :(x[29] == x[13] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:y_idx] == 92
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:id] == 88
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:lifted_constr_ref] == :(x[92] == x[2] * x[3] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:y_idx] == 20
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:id] == 16
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:lifted_constr_ref] == :(x[20] == x[1] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:y_idx] == 109
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:id] == 105
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[109] == x[14] * x[10] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:y_idx] == 93
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:id] == 89
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:lifted_constr_ref] == :(x[93] == x[1] * x[92])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:y_idx] == 85
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:id] == 81
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[85] == x[2] * x[3] * x[4])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:y_idx] == 16
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:id] == 12
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:lifted_constr_ref] == :(x[16] == x[1] * x[11])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:y_idx] == 25
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:id] == 21
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:lifted_constr_ref] == :(x[25] == x[4] * x[21])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:y_idx] == 73
    @test m.internalModel.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:id] == 69
    @test m.internalModel.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:lifted_constr_ref] == :(x[73] == x[53] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:y_idx] == 34
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:id] == 30
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:lifted_constr_ref] == :(x[34] == x[5] * x[32])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:y_idx] == 40
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:id] == 36
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:lifted_constr_ref] == :(x[40] == x[3] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:y_idx] == 67
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:id] == 63
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[67] == x[5] * x[9] * x[10] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:y_idx] == 102
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:id] == 98
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:lifted_constr_ref] == :(x[102] == x[5] * x[9] * x[43])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:y_idx] == 27
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:id] == 23
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:lifted_constr_ref] == :(x[27] == x[4] * x[26])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:y_idx] == 56
    @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:id] == 52
    @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:lifted_constr_ref] == :(x[56] == x[20] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:y_idx] == 63
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:id] == 59
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:lifted_constr_ref] == :(x[63] == x[1] * x[3] * x[4] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:y_idx] == 15
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:id] == 11
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:lifted_constr_ref] == :(x[15] == x[14] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:y_idx] == 77
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:id] == 73
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:lifted_constr_ref] == :(x[77] == x[1] * x[4] * x[76])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:y_idx] == 89
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:id] == 85
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:lifted_constr_ref] == :(x[89] == x[1] * x[88])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:y_idx] == 94
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:id] == 90
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:lifted_constr_ref] == :(x[94] == x[5] * x[90])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:y_idx] == 108
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:id] == 104
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[108] == x[6] * x[10] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:y_idx] == 70
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:id] == 66
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:lifted_constr_ref] == :(x[70] == x[4] * x[48])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:y_idx] == 54
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:id] == 50
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:lifted_constr_ref] == :(x[54] == x[6] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:y_idx] == 99
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:id] == 95
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:lifted_constr_ref] == :(x[99] == x[1] * x[9] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:y_idx] == 36
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:id] == 32
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:lifted_constr_ref] == :(x[36] == x[1] * x[35])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:y_idx] == 14
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:id] == 10
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:lifted_constr_ref] == :(x[14] == x[5] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:y_idx] == 22
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:id] == 18
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:lifted_constr_ref] == :(x[22] == x[4] * x[13])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:y_idx] == 64
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:id] == 60
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[64] == x[1] * x[2] * x[10] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:y_idx] == 72
    @test m.internalModel.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:id] == 68
    @test m.internalModel.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:lifted_constr_ref] == :(x[72] == x[46] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:y_idx] == 55
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:id] == 51
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:lifted_constr_ref] == :(x[55] == x[17] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:y_idx] == 98
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:id] == 94
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:lifted_constr_ref] == :(x[98] == x[2] * x[5] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:y_idx] == 43
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:id] == 39
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[43] == x[10] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:y_idx] == 78
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:id] == 74
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:lifted_constr_ref] == :(x[78] == x[2] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:y_idx] == 52
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:id] == 48
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:lifted_constr_ref] == :(x[52] == x[2] * x[5] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:y_idx] == 91
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:id] == 87
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:lifted_constr_ref] == :(x[91] == x[1] * x[90])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:y_idx] == 104
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:id] == 100
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:lifted_constr_ref] == :(x[104] == x[3] * x[4] * x[17])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:y_idx] == 76
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:id] == 72
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:lifted_constr_ref] == :(x[76] == x[3] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:y_idx] == 79
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:id] == 75
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:lifted_constr_ref] == :(x[79] == x[1] * x[4] * x[78])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:y_idx] == 10
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:id] == 6
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:lifted_constr_ref] == :(x[10] == (*)(x[3]))
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:y_idx] == 71
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:id] == 67
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:lifted_constr_ref] == :(x[71] == x[4] * x[49])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:y_idx] == 88
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:id] == 84
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:lifted_constr_ref] == :(x[88] == x[3] * x[4] * x[9])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:y_idx] == 30
    @test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:id] == 26
    @test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:lifted_constr_ref] == :(x[30] == x[15] * x[28])
    @test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:y_idx] == 68
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:id] == 64
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:lifted_constr_ref] == :(x[68] == x[4] * x[46])
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:y_idx] == 12
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:id] == 8
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:lifted_constr_ref] == :(x[12] == x[5] * x[11])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:y_idx] == 74
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:id] == 70
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:lifted_constr_ref] == :(x[74] == x[1] * x[4] * x[7])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:y_idx] == 105
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:id] == 101
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:lifted_constr_ref] == :(x[105] == x[3] * x[4] * x[20])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:y_idx] == 82
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:id] == 78
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:lifted_constr_ref] == :(x[82] == x[1] * x[4] * x[11])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:y_idx] == 44
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:id] == 40
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:lifted_constr_ref] == :(x[44] == x[9] * x[43])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:y_idx] == 61
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 57
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[61] == x[1] * x[2] * x[3] * x[4])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:y_idx] == 18
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:id] == 14
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:lifted_constr_ref] == :(x[18] == x[17] * x[10])
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:y_idx] == 8
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:id] == 4
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:lifted_constr_ref] == :(x[8] == x[1] * x[7])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 6
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:id] == 2
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:lifted_constr_ref] == :(x[6] == x[1] * x[2])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:y_idx] == 101
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:id] == 97
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:lifted_constr_ref] == :(x[101] == x[1] * x[2] * x[40])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:y_idx] == 60
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:id] == 56
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:lifted_constr_ref] == :(x[60] == x[17] * x[37])
    @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:y_idx] == 103
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:id] == 99
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:lifted_constr_ref] == :(x[103] == x[3] * x[4] * x[6])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:y_idx] == 7
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:id] == 3
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:lifted_constr_ref] == :(x[7] == x[2] * x[3])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:y_idx] == 47
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:id] == 43
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:lifted_constr_ref] == :(x[47] == x[2] * x[3] * x[5])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:nonlinear_type] == :MULTILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:y_idx] == 13
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:id] == 9
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:lifted_constr_ref] == :(x[13] == x[3] * x[6])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:y_idx] == 35
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:id] == 31
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:lifted_constr_ref] == :(x[35] == x[9] * x[31])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:y_idx] == 33
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:id] == 29
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:lifted_constr_ref] == :(x[33] == x[1] * x[32])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:y_idx] == 21
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:id] == 17
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:lifted_constr_ref] == :(x[21] == x[3] * x[20])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:y_idx] == 45
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:id] == 41
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:lifted_constr_ref] == :(x[45] == x[5] * x[44])
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:nonlinear_type] == :BILINEAR

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

@testset "Expression Parsing || corner cases" begin
    @testset "Corner Cases - 1 : sign convertor special case" begin
        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                minlp_solver=pavito_solver,
                                loglevel=100)

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
        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                loglevel=100)

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

    @testset "Corner Cases - 2 : full sub-expression" begin
        test_solver = AlpineSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                loglevel=100)

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

@testset "Expression Parsing || Discrete Multilinear" begin

    @testset "Expression Parsing || bmpl && binlin && binprod" begin

        test_solver=AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)

        m = bpml(solver=test_solver)

        JuMP.build(m) # Setup internal model

        @test length(keys(m.internalModel.nonconvex_terms)) == 12

        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:y_idx] == 11
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 12
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[6])]][:y_idx] == 13
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:y_idx] == 14
        @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[6])]][:y_idx] == 15
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:y_idx] == 16
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[16])]][:y_idx] == 17
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[8])]][:y_idx] == 18
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[18])]][:y_idx] == 19
        @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[18])]][:y_idx] == 20
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]][:y_idx] == 21
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[9]), :(x[10]), :(x[6])]][:y_idx] == 22

        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :BINPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[6])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] == :BINPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[6])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[16])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[8])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[18])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[18])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]][:nonlinear_type] == :BINPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[9]), :(x[10]), :(x[6])]][:nonlinear_type] == :MULTILINEAR

        @test length(m.internalModel.var_type) == 22
        @test m.internalModel.var_type[11] == :Cont
        @test m.internalModel.var_type[12] == :Bin
        @test m.internalModel.var_type[13] == :Cont
        @test m.internalModel.var_type[14] == :Bin
        @test m.internalModel.var_type[15] == :Cont
        @test m.internalModel.var_type[16] == :Cont
        @test m.internalModel.var_type[17] == :Cont
        @test m.internalModel.var_type[18] == :Cont
        @test m.internalModel.var_type[19] == :Cont
        @test m.internalModel.var_type[20] == :Cont
        @test m.internalModel.var_type[21] == :Bin
        @test m.internalModel.var_type[22] == :Cont
    end

    @testset "Expression Parsing || bmpl && binlin && binprod with linear lifting and coefficients" begin

        test_solver=AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)

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
        nlk8 = Expr[:(x[24]), :(x[23])]
        nlk9 = Expr[:(x[1]), :(x[2])]
        nlk10 = Expr[:(x[2]), :(x[3])]
        nlk11 = Expr[:(x[20]), :(x[21])]
        nlk12 = Expr[:(x[3]), :(x[4])]

        @test haskey(m.internalModel.nonconvex_terms, nlk1)
        @test haskey(m.internalModel.nonconvex_terms, nlk2)
        @test haskey(m.internalModel.nonconvex_terms, nlk3)
        @test haskey(m.internalModel.nonconvex_terms, nlk4)
        @test haskey(m.internalModel.nonconvex_terms, nlk5)
        @test haskey(m.internalModel.nonconvex_terms, nlk6)
        @test haskey(m.internalModel.nonconvex_terms, nlk7)
        @test haskey(m.internalModel.nonconvex_terms, nlk8)
        @test haskey(m.internalModel.nonconvex_terms, nlk9)
        @test haskey(m.internalModel.nonconvex_terms, nlk10)
        @test haskey(m.internalModel.nonconvex_terms, nlk11)
        @test haskey(m.internalModel.nonconvex_terms, nlk12)

        @test m.internalModel.nonconvex_terms[nlk1][:id] == 7
        @test m.internalModel.nonconvex_terms[nlk2][:id] == 1
        @test m.internalModel.nonconvex_terms[nlk3][:id] == 9
        @test m.internalModel.nonconvex_terms[nlk4][:id] == 12
        @test m.internalModel.nonconvex_terms[nlk5][:id] == 3
        @test m.internalModel.nonconvex_terms[nlk6][:id] == 4
        @test m.internalModel.nonconvex_terms[nlk7][:id] == 5
        @test m.internalModel.nonconvex_terms[nlk8][:id] == 11
        @test m.internalModel.nonconvex_terms[nlk9][:id] == 6
        @test m.internalModel.nonconvex_terms[nlk10][:id] == 2
        @test m.internalModel.nonconvex_terms[nlk11][:id] == 8
        @test m.internalModel.nonconvex_terms[nlk12][:id] == 10

        @test m.internalModel.nonconvex_terms[nlk1][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk2][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk3][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk4][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk5][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk6][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk7][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk8][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk9][:y_type] == :Bin
        @test m.internalModel.nonconvex_terms[nlk10][:y_type] == :Bin
        @test m.internalModel.nonconvex_terms[nlk11][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[nlk12][:y_type] == :Bin

        @test m.internalModel.nonconvex_terms[nlk1][:y_idx] == 19
        @test m.internalModel.nonconvex_terms[nlk2][:y_idx] == 11
        @test m.internalModel.nonconvex_terms[nlk3][:y_idx] == 23
        @test m.internalModel.nonconvex_terms[nlk4][:y_idx] == 27
        @test m.internalModel.nonconvex_terms[nlk5][:y_idx] == 14
        @test m.internalModel.nonconvex_terms[nlk6][:y_idx] == 16
        @test m.internalModel.nonconvex_terms[nlk7][:y_idx] == 17
        @test m.internalModel.nonconvex_terms[nlk8][:y_idx] == 25
        @test m.internalModel.nonconvex_terms[nlk9][:y_idx] == 18
        @test m.internalModel.nonconvex_terms[nlk10][:y_idx] == 12
        @test m.internalModel.nonconvex_terms[nlk11][:y_idx] == 22
        @test m.internalModel.nonconvex_terms[nlk12][:y_idx] == 24

        @test m.internalModel.nonconvex_terms[nlk1][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[nlk2][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[nlk3][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk4][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[nlk5][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[nlk6][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[nlk7][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[nlk8][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[nlk9][:nonlinear_type] == :BINPROD
        @test m.internalModel.nonconvex_terms[nlk10][:nonlinear_type] == :BINPROD
        @test m.internalModel.nonconvex_terms[nlk11][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[nlk12][:nonlinear_type] == :BINPROD
    end

    @testset "Expression Parsing || INTPROD Operators" begin
        test_solver=AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0),mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = intprod_basic(solver=test_solver)

        JuMP.build(m) # Setup internal model

        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:y_idx] == 18
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:id] == 4
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[18] == x[2] * x[3] * x[4])
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:constr_id] == Set(Any[3])
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[16])]][:y_idx] == 17
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[16])]][:id] == 3
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[16])]][:lifted_constr_ref] == :(x[17] == x[1] * x[4] * x[16])
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[16])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[16])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[16])]][:constr_id] == Set(Any[2])
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[2])]][:y_idx] == 12
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[2])]][:id] == 1
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[2])]][:lifted_constr_ref] == :(x[12] == x[11] * x[2])
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[2])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[2])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[2])]][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:y_idx] == 19
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:id] == 5
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:lifted_constr_ref] == :(x[19] == (*)(x[5]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:constr_id] == Set(Any[4])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[20]), :(x[19])]][:y_idx] == 21
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[20]), :(x[19])]][:id] == 6
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[20]), :(x[19])]][:lifted_constr_ref] == :(x[21] == x[5] * x[20] * x[19])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[20]), :(x[19])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[20]), :(x[19])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[20]), :(x[19])]][:constr_id] == Set(Any[4])
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:y_idx] == 15
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:id] == 2
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:lifted_constr_ref] == :(x[15] == x[13] * x[14])
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:constr_id] == Set(Any[1])
        lk1 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 4), (1.0, 1), (1.0, 5), (1.0, 3)])))
        @test m.internalModel.linear_terms[lk1][:y_idx] == 11
        @test m.internalModel.linear_terms[lk1][:id] == 1
        @test m.internalModel.linear_terms[lk1][:y_type] == :(Int)
        lk2 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (-1.0, 3)])))
        @test m.internalModel.linear_terms[lk2][:y_idx] == 20
        @test m.internalModel.linear_terms[lk2][:id] == 5
        @test m.internalModel.linear_terms[lk2][:y_type] == :(Int)
        lk3 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 3)])))
        @test m.internalModel.linear_terms[lk3][:y_idx] == 16
        @test m.internalModel.linear_terms[lk3][:id] == 4
        @test m.internalModel.linear_terms[lk3][:y_type] == :(Int)
        lk4 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4), (1.0, 3)])))
        @test m.internalModel.linear_terms[lk4][:y_idx] == 14
        @test m.internalModel.linear_terms[lk4][:id] == 3
        @test m.internalModel.linear_terms[lk4][:y_type] == :(Int)
        lk5 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1)])))
        @test m.internalModel.linear_terms[lk5][:y_idx] == 13
        @test m.internalModel.linear_terms[lk5][:id] == 2
        @test m.internalModel.linear_terms[lk5][:y_type] == :(Int)
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[15])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[2][:rhs] == 25.0
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[17])]
        @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[3][:rhs] == 10.0
        @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[1]), :(x[18])]
        @test m.internalModel.bounding_constr_mip[3][:coefs] == Any[1.0, -1.0]
        @test m.internalModel.bounding_constr_mip[3][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[3][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[4][:rhs] == 40.0
        @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[21])]
        @test m.internalModel.bounding_constr_mip[4][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[4][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[4][:cnt] == 1
    end

    @testset "Expression Parsing || ex1225a" begin
        test_solver = AlpineSolver(minlp_solver=pavito_solver, nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = ex1225a(solver=test_solver)

        JuMP.build(m)

        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44])]][:y_idx] == 46
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44])]][:id] == 16
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44])]][:lifted_constr_ref] == :(x[46] == (*)(x[44]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44])]][:constr_id] == Set(Any[25, 28])
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[52])]][:y_idx] == 53
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[52])]][:id] == 22
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[52])]][:lifted_constr_ref] == :(x[53] == x[13] * x[52])
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[52])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[52])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[52])]][:constr_id] == Set(Any[26])
        @test m.internalModel.nonconvex_terms[Expr[:(x[25]), :(x[10])]][:y_idx] == 64
        @test m.internalModel.nonconvex_terms[Expr[:(x[25]), :(x[10])]][:id] == 33
        @test m.internalModel.nonconvex_terms[Expr[:(x[25]), :(x[10])]][:lifted_constr_ref] == :(x[64] == x[25] * x[10])
        @test m.internalModel.nonconvex_terms[Expr[:(x[25]), :(x[10])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[25]), :(x[10])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[25]), :(x[10])]][:constr_id] == Set(Any[35])
        @test m.internalModel.nonconvex_terms[Expr[:(x[19]), :(x[36])]][:y_idx] == 37
        @test m.internalModel.nonconvex_terms[Expr[:(x[19]), :(x[36])]][:id] == 9
        @test m.internalModel.nonconvex_terms[Expr[:(x[19]), :(x[36])]][:lifted_constr_ref] == :(x[37] == x[19] * x[36])
        @test m.internalModel.nonconvex_terms[Expr[:(x[19]), :(x[36])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[19]), :(x[36])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[19]), :(x[36])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38])]][:y_idx] == 40
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38])]][:id] == 11
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38])]][:lifted_constr_ref] == :(x[40] == (*)(x[38]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38])]][:constr_id] == Set(Any[27, 24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[35]), :(x[34])]][:y_idx] == 36
        @test m.internalModel.nonconvex_terms[Expr[:(x[35]), :(x[34])]][:id] == 8
        @test m.internalModel.nonconvex_terms[Expr[:(x[35]), :(x[34])]][:lifted_constr_ref] == :(x[36] == x[35] * x[34])
        @test m.internalModel.nonconvex_terms[Expr[:(x[35]), :(x[34])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[35]), :(x[34])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[35]), :(x[34])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[23])]][:y_idx] == 27
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[23])]][:id] == 1
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[23])]][:lifted_constr_ref] == :(x[27] == x[20] * x[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[23])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[23])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[23])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[48])]][:y_idx] == 49
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[48])]][:id] == 19
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[48])]][:lifted_constr_ref] == :(x[49] == x[6] * x[48])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[48])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[48])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[48])]][:constr_id] == Set(Any[25])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:y_idx] == 56
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:id] == 25
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:lifted_constr_ref] == :(x[56] == x[5] * x[11])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:constr_id] == Set(Any[27])
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[11])]][:y_idx] == 42
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[11])]][:id] == 13
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[11])]][:lifted_constr_ref] == :(x[42] == (*)(x[11]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[11])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[11])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[11])]][:constr_id] == Set(Any[27, 24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38]), :(x[38])]][:y_idx] == 39
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38]), :(x[38])]][:id] == 10
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38]), :(x[38])]][:lifted_constr_ref] == :(x[39] == (*)(x[38]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38]), :(x[38])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38]), :(x[38])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[38]), :(x[38]), :(x[38])]][:constr_id] == Set(Any[24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[13])]][:y_idx] == 58
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[13])]][:id] == 27
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[13])]][:lifted_constr_ref] == :(x[58] == x[7] * x[13])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[13])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[13])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[13])]][:constr_id] == Set(Any[29])
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[12])]][:y_idx] == 48
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[12])]][:id] == 18
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[12])]][:lifted_constr_ref] == :(x[48] == (*)(x[12]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[12])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[12])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[12])]][:constr_id] == Set(Any[25, 28])
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[11])]][:y_idx] == 59
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[11])]][:id] == 28
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[11])]][:lifted_constr_ref] == :(x[59] == x[20] * x[11])
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[11])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[11])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[20]), :(x[11])]][:constr_id] == Set(Any[30])
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[40])]][:y_idx] == 41
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[40])]][:id] == 12
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[40])]][:lifted_constr_ref] == :(x[41] == x[11] * x[40])
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[40])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[40])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[40])]][:constr_id] == Set(Any[24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[25])]][:y_idx] == 35
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[25])]][:id] == 7
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[25])]][:lifted_constr_ref] == :(x[35] == x[22] * x[25])
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[25])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[25])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[25])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[28])]][:y_idx] == 29
        @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[28])]][:id] == 3
        @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[28])]][:lifted_constr_ref] == :(x[29] == x[17] * x[28])
        @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[28])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[28])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[17]), :(x[28])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[18]), :(x[32])]][:y_idx] == 33
        @test m.internalModel.nonconvex_terms[Expr[:(x[18]), :(x[32])]][:id] == 6
        @test m.internalModel.nonconvex_terms[Expr[:(x[18]), :(x[32])]][:lifted_constr_ref] == :(x[33] == x[18] * x[32])
        @test m.internalModel.nonconvex_terms[Expr[:(x[18]), :(x[32])]][:nonlinear_type] == :BINLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[18]), :(x[32])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[18]), :(x[32])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50])]][:y_idx] == 52
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50])]][:id] == 21
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50])]][:lifted_constr_ref] == :(x[52] == (*)(x[50]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50])]][:constr_id] == Set(Any[26, 29])
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[13])]][:y_idx] == 61
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[13])]][:id] == 30
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[13])]][:lifted_constr_ref] == :(x[61] == x[22] * x[13])
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[13])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[13])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[22]), :(x[13])]][:constr_id] == Set(Any[32])
        @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[26])]][:y_idx] == 28
        @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[26])]][:id] == 2
        @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[26])]][:lifted_constr_ref] == :(x[28] == x[27] * x[26])
        @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[26])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[26])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[26])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[54])]][:y_idx] == 55
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[54])]][:id] == 24
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[54])]][:lifted_constr_ref] == :(x[55] == x[7] * x[54])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[54])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[54])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[54])]][:constr_id] == Set(Any[26])
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44]), :(x[44])]][:y_idx] == 45
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44]), :(x[44])]][:id] == 15
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44]), :(x[44])]][:lifted_constr_ref] == :(x[45] == (*)(x[44]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44]), :(x[44])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44]), :(x[44])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[44]), :(x[44])]][:constr_id] == Set(Any[25])
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[24])]][:y_idx] == 31
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[24])]][:id] == 4
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[24])]][:lifted_constr_ref] == :(x[31] == x[21] * x[24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[24])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[24])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[24])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[8])]][:y_idx] == 62
        @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[8])]][:id] == 31
        @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[8])]][:lifted_constr_ref] == :(x[62] == x[23] * x[8])
        @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[8])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[8])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[8])]][:constr_id] == Set(Any[33])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[12])]][:y_idx] == 57
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[12])]][:id] == 26
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[12])]][:lifted_constr_ref] == :(x[57] == x[6] * x[12])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[12])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[12])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[12])]][:constr_id] == Set(Any[28])
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[46])]][:y_idx] == 47
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[46])]][:id] == 17
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[46])]][:lifted_constr_ref] == :(x[47] == x[12] * x[46])
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[46])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[46])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[46])]][:constr_id] == Set(Any[25])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[42])]][:y_idx] == 43
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[42])]][:id] == 14
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[42])]][:lifted_constr_ref] == :(x[43] == x[5] * x[42])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[42])]][:nonlinear_type] == :BILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[42])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[42])]][:constr_id] == Set(Any[24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[13])]][:y_idx] == 54
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[13])]][:id] == 23
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[13])]][:lifted_constr_ref] == :(x[54] == (*)(x[13]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[13])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[13])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[13])]][:constr_id] == Set(Any[26, 29])
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[12])]][:y_idx] == 60
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[12])]][:id] == 29
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[12])]][:lifted_constr_ref] == :(x[60] == x[21] * x[12])
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[12])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[12])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[21]), :(x[12])]][:constr_id] == Set(Any[31])
        @test m.internalModel.nonconvex_terms[Expr[:(x[24]), :(x[9])]][:y_idx] == 63
        @test m.internalModel.nonconvex_terms[Expr[:(x[24]), :(x[9])]][:id] == 32
        @test m.internalModel.nonconvex_terms[Expr[:(x[24]), :(x[9])]][:lifted_constr_ref] == :(x[63] == x[24] * x[9])
        @test m.internalModel.nonconvex_terms[Expr[:(x[24]), :(x[9])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[24]), :(x[9])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[24]), :(x[9])]][:constr_id] == Set(Any[34])
        @test m.internalModel.nonconvex_terms[Expr[:(x[31]), :(x[30])]][:y_idx] == 32
        @test m.internalModel.nonconvex_terms[Expr[:(x[31]), :(x[30])]][:id] == 5
        @test m.internalModel.nonconvex_terms[Expr[:(x[31]), :(x[30])]][:lifted_constr_ref] == :(x[32] == x[31] * x[30])
        @test m.internalModel.nonconvex_terms[Expr[:(x[31]), :(x[30])]][:nonlinear_type] == :INTLIN
        @test m.internalModel.nonconvex_terms[Expr[:(x[31]), :(x[30])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[31]), :(x[30])]][:constr_id] == Set(Any[23])
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50]), :(x[50])]][:y_idx] == 51
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50]), :(x[50])]][:id] == 20
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50]), :(x[50])]][:lifted_constr_ref] == :(x[51] == (*)(x[50]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50]), :(x[50])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50]), :(x[50])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[50]), :(x[50]), :(x[50])]][:constr_id] == Set(Any[26])
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 1.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[14]), :(x[15]), :(x[16])]
        @test isapprox(m.internalModel.bounding_constr_mip[1][:coefs],Any[1.0, 1.0, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 3
        @test m.internalModel.bounding_constr_mip[2][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[5]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[2][:coefs],Any[0.000338983, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[3][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[6]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[3][:coefs],Any[0.000338983, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[3][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[3][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[4][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[7]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[4][:coefs],Any[0.000338983, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[4][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[4][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[5][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[5][:vars] == Any[:(x[2]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[5][:coefs],Any[0.0125, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[5][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[5][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[6][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[6][:vars] == Any[:(x[3]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[6][:coefs],Any[0.04, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[6][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[6][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[7][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[7][:vars] == Any[:(x[4]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[7][:coefs],Any[0.0222222, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[7][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[7][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[8][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[8][:vars] == Any[:(x[8]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[8][:coefs],Any[0.0025, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[8][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[8][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[9][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[9][:vars] == Any[:(x[9]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[9][:coefs],Any[0.0025, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[9][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[9][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[10][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[10][:vars] == Any[:(x[10]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[10][:coefs],Any[0.0025, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[10][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[10][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[11][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[11][:vars] == Any[:(x[11]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[11][:coefs],Any[0.00285714, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[11][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[11][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[12][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[12][:vars] == Any[:(x[12]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[12][:coefs],Any[0.00285714, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[12][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[12][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[13][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[13][:vars] == Any[:(x[13]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[13][:coefs],Any[0.00285714, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[13][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[13][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[14][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[14][:vars] == Any[:(x[14]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[14][:coefs],Any[1.0, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[14][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[14][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[15][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[15][:vars] == Any[:(x[15]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[15][:coefs],Any[1.0, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[15][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[15][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[16][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[16][:vars] == Any[:(x[16]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[16][:coefs],Any[1.0, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[16][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[16][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[17][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[17][:vars] == Any[:(x[20]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[17][:coefs],Any[1.0, -3.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[17][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[17][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[18][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[18][:vars] == Any[:(x[21]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[18][:coefs],Any[1.0, -3.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[18][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[18][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[19][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[19][:vars] == Any[:(x[22]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[19][:coefs],Any[1.0, -3.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[19][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[19][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[20][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[20][:vars] == Any[:(x[23]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[20][:coefs],Any[1.0, -3.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[20][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[20][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[21][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[21][:vars] == Any[:(x[24]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[21][:coefs],Any[1.0, -3.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[21][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[21][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[22][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[22][:vars] == Any[:(x[25]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[22][:coefs],Any[1.0, -3.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[22][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[22][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[23][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[23][:vars] == Any[:(x[29]), :(x[33]), :(x[37]), :(x[1])]
        @test isapprox(m.internalModel.bounding_constr_mip[23][:coefs],Any[-1.0, -1.0, -1.0, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[23][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[23][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[24][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[24][:vars] == Any[:(x[39]), :(x[41]), :(x[43]), :(x[2])]
        @test isapprox(m.internalModel.bounding_constr_mip[24][:coefs],Any[-19.9, -0.161, 1.90169e-7, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[24][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[24][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[25][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[25][:vars] == Any[:(x[45]), :(x[47]), :(x[49]), :(x[3])]
        @test isapprox(m.internalModel.bounding_constr_mip[25][:coefs],Any[-1.21, -0.0644, 1.91186e-7, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[25][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[25][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[26][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[26][:vars] == Any[:(x[51]), :(x[53]), :(x[55]), :(x[4])]
        @test isapprox(m.internalModel.bounding_constr_mip[26][:coefs],Any[-6.52, -0.102, 7.86441e-8, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[26][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[26][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[27][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[27][:vars] == Any[:(x[56]), :(x[40]), :(x[42]), :(x[8])]
        @test isapprox(m.internalModel.bounding_constr_mip[27][:coefs],Any[-0.000235932, -629.0, 0.0116, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[27][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[27][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[28][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[28][:vars] == Any[:(x[57]), :(x[46]), :(x[48]), :(x[9])]
        @test isapprox(m.internalModel.bounding_constr_mip[28][:coefs],Any[-0.001, -215.0, 0.115, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[28][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[28][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[29][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[29][:vars] == Any[:(x[58]), :(x[52]), :(x[54]), :(x[10])]
        @test isapprox(m.internalModel.bounding_constr_mip[29][:coefs],Any[-0.000179661, -361.0, 0.00946, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[29][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[29][:cnt] == 4
        @test m.internalModel.bounding_constr_mip[30][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[30][:vars] == Any[:(x[59]), :(x[14])]
        @test isapprox(m.internalModel.bounding_constr_mip[30][:coefs],Any[0.00285714, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[30][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[30][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[31][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[31][:vars] == Any[:(x[60]), :(x[15])]
        @test isapprox(m.internalModel.bounding_constr_mip[31][:coefs],Any[0.00285714, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[31][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[31][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[32][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[32][:vars] == Any[:(x[61]), :(x[16])]
        @test isapprox(m.internalModel.bounding_constr_mip[32][:coefs],Any[0.00285714, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[32][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[32][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[33][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[33][:vars] == Any[:(x[62]), :(x[17])]
        @test isapprox(m.internalModel.bounding_constr_mip[33][:coefs],Any[0.0025, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[33][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[33][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[34][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[34][:vars] == Any[:(x[63]), :(x[18])]
        @test isapprox(m.internalModel.bounding_constr_mip[34][:coefs],Any[0.0025, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[34][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[34][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[35][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[35][:vars] == Any[:(x[64]), :(x[19])]
        @test isapprox(m.internalModel.bounding_constr_mip[35][:coefs],Any[0.0025, -1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[35][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[35][:cnt] == 2
    end

    @testset "Expression Parsing || prob10" begin
        test_solver = AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = prob10(solver=test_solver)

        JuMP.build(m)

        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:y_idx] == 7
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:id] == 2
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:lifted_constr_ref] == :(x[7] == (*)(x[6]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:constr_id] == Set(Any[3])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 9
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 3
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[9] == sin(x[8]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[3])
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:y_idx] == 5
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:id] == 1
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:lifted_constr_ref] == :(x[5] == (*)(x[4]))
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:constr_id] == Set(Any[3])
        lk1 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7), (1.0, 5)])))
        @test m.internalModel.linear_terms[lk1][:y_idx] == 8
        @test m.internalModel.linear_terms[lk1][:id] == 3
        @test m.internalModel.linear_terms[lk1][:y_type] == :(Cont)
        lk2 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, -10.0),Pair{Symbol,Any}(:coef_var, Set(Any[(2.0, 2)])))
        @test m.internalModel.linear_terms[lk2][:y_idx] == 4
        @test m.internalModel.linear_terms[lk2][:id] == 1
        @test m.internalModel.linear_terms[lk2][:y_type] == :(Cont)
        lk3 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, -5.0),Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3)])))
        @test m.internalModel.linear_terms[lk3][:y_idx] == 6
        @test m.internalModel.linear_terms[lk3][:id] == 2
        @test m.internalModel.linear_terms[lk3][:y_type] == :(Int)
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 7.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[2]), :(x[3])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[0.7, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[2][:rhs] == 19.0
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[2]), :(x[3])]
        @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[2.5, 1.0]
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 2
        @test m.internalModel.bounding_constr_mip[3][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[5]), :(x[7]), :(x[9]), :(x[1])]
        @test m.internalModel.bounding_constr_mip[3][:coefs] == Any[1.1, 1.1, 1.0, -1.0]
        @test m.internalModel.bounding_constr_mip[3][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[3][:cnt] == 4
    end

    @testset "Expression Parsing || prob03" begin
        test_solver = AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = prob03(solver=test_solver)

        JuMP.build(m)

        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:y_idx] == 4
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:id] == 1
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:lifted_constr_ref] == :(x[4] == x[2] * x[3])
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:nonlinear_type] == :INTPROD
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:y_type] == :Int
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:constr_id] == Set(Any[2])
        @test m.internalModel.bounding_constr_mip[1][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[2]), :(x[3]), :(x[1])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[-3.0, -2.0, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 3
        @test m.internalModel.bounding_constr_mip[2][:rhs] == -3.5
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[4])]
        @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[-1.0]
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 1
    end

    @testset "Expression Parsing || st_miqp5" begin
        test_solver = AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = st_miqp5(solver=test_solver)

        JuMP.build(m)

        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:y_idx] == 10
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:id] == 2
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:lifted_constr_ref] == :(x[10] == x[6] * x[6])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[6])]][:constr_id] == Set(Any[14])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:y_idx] == 9
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:id] == 1
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:lifted_constr_ref] == :(x[9] == x[5] * x[5])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:nonlinear_type] == :MONOMIAL
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[5])]][:constr_id] == Set(Any[14])
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 60.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[1][:coefs],Any[-1.93415, 1.80315, 2.89696, 0.729325, 3.88374];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[2][:rhs] == 60.0
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[2][:coefs],Any[-1.13151, 1.10501, -1.01839, 2.62557, 4.85468];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[3][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[3][:coefs],Any[-0.05248, -0.904838, 0.209521, -0.29173, -0.222506];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[3][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[3][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[4][:rhs] == 1.0
        @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[4][:coefs],Any[0.05248, 0.904838, -0.209521, 0.29173, 0.222506];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[4][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[4][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[5][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[5][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[5][:coefs],Any[0.445392, 0.30152, 0.587645, -0.145865, -0.586607];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[5][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[5][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[6][:rhs] == 1.0
        @test m.internalModel.bounding_constr_mip[6][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[6][:coefs],Any[-0.445392, -0.30152, -0.587645, 0.145865, 0.586607];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[6][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[6][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[7][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[7][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[7][:coefs],Any[-0.328189, 0.199987, 0.506106, -0.58346, 0.505696];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[7][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[7][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[8][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[8][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[8][:coefs],Any[-0.345682, -0.101626, 0.575947, 0.729325, 0.0809113];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[8][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[8][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[9][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[9][:vars] == Any[:(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[9][:coefs],Any[0.756087, -0.200079, 0.151379, 0.145865, 0.586607];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[9][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[9][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[10][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[10][:vars] == Any[:(x[7]), :(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[10][:coefs],Any[-1.0, 0.05248, 0.904838, -0.209521, 0.29173, 0.222506];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[10][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[10][:cnt] == 6
        @test m.internalModel.bounding_constr_mip[11][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[11][:vars] == Any[:(x[7]), :(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[11][:coefs],Any[1.0, -0.05248, -0.904838, 0.209521, -0.29173, -0.222506];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[11][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[11][:cnt] == 6
        @test m.internalModel.bounding_constr_mip[12][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[12][:vars] == Any[:(x[8]), :(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[12][:coefs],Any[-1.0, -0.445392, -0.30152, -0.587645, 0.145865, 0.586607];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[12][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[12][:cnt] == 6
        @test m.internalModel.bounding_constr_mip[13][:rhs] == -0.0
        @test m.internalModel.bounding_constr_mip[13][:vars] == Any[:(x[8]), :(x[2]), :(x[3]), :(x[4]), :(x[5]), :(x[6])]
        @test isapprox(m.internalModel.bounding_constr_mip[13][:coefs],Any[1.0, 0.445392, 0.30152, 0.587645, -0.145865, -0.586607];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[13][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[13][:cnt] == 6
        @test m.internalModel.bounding_constr_mip[14][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[14][:vars] == Any[:(x[9]), :(x[5]), :(x[10]), :(x[6]), :(x[2]), :(x[3]), :(x[4]), :(x[1])]
        @test isapprox(m.internalModel.bounding_constr_mip[14][:coefs],Any[-5.0, 0.87519, -52.0, 192.711, 54.0616, 45.2691, 33.0896, 1.0];atol=1e-3)
        @test m.internalModel.bounding_constr_mip[14][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[14][:cnt] == 8
    end

    @testset "Expression Parsing || discretemulti_basic" begin

		test_solver = AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
		m = discretemulti_basic(solver=test_solver)

		JuMP.build(m)

		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:y_idx] == 18
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:id] == 3
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:lifted_constr_ref] == :(x[18] == x[3] * x[8])
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:constr_id] == Set(Any[0])
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[9])]][:y_idx] == 33
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[9])]][:id] == 18
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[9])]][:lifted_constr_ref] == :(x[33] == x[32] * x[9])
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[9])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[9])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[9])]][:constr_id] == Set(Any[12])
		@test m.internalModel.nonconvex_terms[Expr[:(x[40]), :(x[41])]][:y_idx] == 42
		@test m.internalModel.nonconvex_terms[Expr[:(x[40]), :(x[41])]][:id] == 27
		@test m.internalModel.nonconvex_terms[Expr[:(x[40]), :(x[41])]][:lifted_constr_ref] == :(x[42] == x[40] * x[41])
		@test m.internalModel.nonconvex_terms[Expr[:(x[40]), :(x[41])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[40]), :(x[41])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[40]), :(x[41])]][:constr_id] == Set(Any[18])
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[8])]][:y_idx] == 23
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[8])]][:id] == 8
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[8])]][:lifted_constr_ref] == :(x[23] == x[13] * x[8])
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[8])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[8])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[8])]][:constr_id] == Set(Any[15, 6])
		@test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[45])]][:y_idx] == 46
		@test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[45])]][:id] == 31
		@test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[45])]][:lifted_constr_ref] == :(x[46] == x[44] * x[45])
		@test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[45])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[45])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[44]), :(x[45])]][:constr_id] == Set(Any[19])
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:y_idx] == 48
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:id] == 33
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:lifted_constr_ref] == :(x[48] == x[3] * x[4])
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:nonlinear_type] == :BINPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:y_type] == :Bin
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:constr_id] == Set(Any[20])
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[9])]][:y_idx] == 24
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[9])]][:id] == 9
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[9])]][:lifted_constr_ref] == :(x[24] == x[14] * x[9])
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[9])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[9])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[9])]][:constr_id] == Set(Any[7, 16])
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:y_idx] == 19
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:id] == 4
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:lifted_constr_ref] == :(x[19] == x[4] * x[9])
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:constr_id] == Set(Any[0])
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[24])]][:y_idx] == 37
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[24])]][:id] == 22
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[24])]][:lifted_constr_ref] == :(x[37] == x[4] * x[24])
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[24])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[24])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[24])]][:constr_id] == Set(Any[16])
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:y_idx] == 16
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:id] == 1
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:lifted_constr_ref] == :(x[16] == x[1] * x[6])
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:constr_id] == Set(Any[0])
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[22])]][:y_idx] == 35
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[22])]][:id] == 20
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[22])]][:lifted_constr_ref] == :(x[35] == x[2] * x[22])
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[22])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[22])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[22])]][:constr_id] == Set(Any[14])
		@test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[10])]][:y_idx] == 25
		@test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[10])]][:id] == 10
		@test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[10])]][:lifted_constr_ref] == :(x[25] == x[15] * x[10])
		@test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[10])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[10])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[15]), :(x[10])]][:constr_id] == Set(Any[17, 8])
		@test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:y_idx] == 51
		@test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:id] == 36
		@test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:lifted_constr_ref] == :(x[51] == x[9] * x[10])
		@test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :BILINEAR
		@test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:constr_id] == Set(Any[21])
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[7])]][:y_idx] == 29
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[7])]][:id] == 14
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[7])]][:lifted_constr_ref] == :(x[29] == x[28] * x[7])
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[7])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[7])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[7])]][:constr_id] == Set(Any[10])
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[8])]][:y_idx] == 31
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[8])]][:id] == 16
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[8])]][:lifted_constr_ref] == :(x[31] == x[30] * x[8])
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[8])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[8])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[8])]][:constr_id] == Set(Any[11])
		@test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:y_idx] == 47
		@test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:id] == 32
		@test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:lifted_constr_ref] == :(x[47] == x[8] * x[9])
		@test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:nonlinear_type] == :BILINEAR
		@test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:constr_id] == Set(Any[20])
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[15])]][:y_idx] == 32
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[15])]][:id] == 17
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[15])]][:lifted_constr_ref] == :(x[32] == x[14] * x[15])
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[15])]][:nonlinear_type] == :INTPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[15])]][:y_type] == :Int
		@test m.internalModel.nonconvex_terms[Expr[:(x[14]), :(x[15])]][:constr_id] == Set(Any[21, 12])
		@test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:y_idx] == 43
		@test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:id] == 28
		@test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:lifted_constr_ref] == :(x[43] == x[7] * x[8])
		@test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:nonlinear_type] == :BILINEAR
		@test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:constr_id] == Set(Any[19])
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[23])]][:y_idx] == 36
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[23])]][:id] == 21
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[23])]][:lifted_constr_ref] == :(x[36] == x[3] * x[23])
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[23])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[23])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[23])]][:constr_id] == Set(Any[15])
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:y_idx] == 17
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:id] == 2
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:lifted_constr_ref] == :(x[17] == x[2] * x[7])
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:constr_id] == Set(Any[0])
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[7])]][:y_idx] == 22
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[7])]][:id] == 7
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[7])]][:lifted_constr_ref] == :(x[22] == x[12] * x[7])
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[7])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[7])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[7])]][:constr_id] == Set(Any[14, 5])
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[12])]][:y_idx] == 26
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[12])]][:id] == 11
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[12])]][:lifted_constr_ref] == :(x[26] == x[11] * x[12])
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[12])]][:nonlinear_type] == :INTPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[12])]][:y_type] == :Int
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[12])]][:constr_id] == Set(Any[9, 18])
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5])]][:y_idx] == 52
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5])]][:id] == 37
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5])]][:lifted_constr_ref] == :(x[52] == x[4] * x[5])
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5])]][:nonlinear_type] == :BINPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5])]][:y_type] == :Bin
		@test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5])]][:constr_id] == Set(Any[21])
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[21])]][:y_idx] == 34
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[21])]][:id] == 19
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[21])]][:lifted_constr_ref] == :(x[34] == x[1] * x[21])
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[21])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[21])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[21])]][:constr_id] == Set(Any[13])
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 40
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:id] == 25
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:lifted_constr_ref] == :(x[40] == x[1] * x[2])
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :BINPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:y_type] == :Bin
		@test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:constr_id] == Set(Any[18])
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:y_idx] == 30
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:id] == 15
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:lifted_constr_ref] == :(x[30] == x[13] * x[14])
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:nonlinear_type] == :INTPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:y_type] == :Int
		@test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[14])]][:constr_id] == Set(Any[11, 20])
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[43])]][:y_idx] == 45
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[43])]][:id] == 30
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[43])]][:lifted_constr_ref] == :(x[45] == x[28] * x[43])
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[43])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[43])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[28]), :(x[43])]][:constr_id] == Set(Any[19])
		@test m.internalModel.nonconvex_terms[Expr[:(x[52]), :(x[53])]][:y_idx] == 54
		@test m.internalModel.nonconvex_terms[Expr[:(x[52]), :(x[53])]][:id] == 39
		@test m.internalModel.nonconvex_terms[Expr[:(x[52]), :(x[53])]][:lifted_constr_ref] == :(x[54] == x[52] * x[53])
		@test m.internalModel.nonconvex_terms[Expr[:(x[52]), :(x[53])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[52]), :(x[53])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[52]), :(x[53])]][:constr_id] == Set(Any[21])
		@test m.internalModel.nonconvex_terms[Expr[:(x[48]), :(x[49])]][:y_idx] == 50
		@test m.internalModel.nonconvex_terms[Expr[:(x[48]), :(x[49])]][:id] == 35
		@test m.internalModel.nonconvex_terms[Expr[:(x[48]), :(x[49])]][:lifted_constr_ref] == :(x[50] == x[48] * x[49])
		@test m.internalModel.nonconvex_terms[Expr[:(x[48]), :(x[49])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[48]), :(x[49])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[48]), :(x[49])]][:constr_id] == Set(Any[20])
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[51])]][:y_idx] == 53
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[51])]][:id] == 38
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[51])]][:lifted_constr_ref] == :(x[53] == x[32] * x[51])
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[51])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[51])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[32]), :(x[51])]][:constr_id] == Set(Any[21])
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:y_idx] == 44
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:id] == 29
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:lifted_constr_ref] == :(x[44] == x[2] * x[3])
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:nonlinear_type] == :BINPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:y_type] == :Bin
		@test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:constr_id] == Set(Any[19])
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[39])]][:y_idx] == 41
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[39])]][:id] == 26
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[39])]][:lifted_constr_ref] == :(x[41] == x[26] * x[39])
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[39])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[39])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[39])]][:constr_id] == Set(Any[18])
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[6])]][:y_idx] == 27
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[6])]][:id] == 12
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[6])]][:lifted_constr_ref] == :(x[27] == x[26] * x[6])
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[6])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[6])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[26]), :(x[6])]][:constr_id] == Set(Any[9])
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:y_idx] == 28
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:id] == 13
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:lifted_constr_ref] == :(x[28] == x[12] * x[13])
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:nonlinear_type] == :INTPROD
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:y_type] == :Int
		@test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:constr_id] == Set(Any[10, 19])
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[25])]][:y_idx] == 38
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[25])]][:id] == 23
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[25])]][:lifted_constr_ref] == :(x[38] == x[5] * x[25])
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[25])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[25])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[25])]][:constr_id] == Set(Any[17])
		@test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:y_idx] == 39
		@test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:id] == 24
		@test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:lifted_constr_ref] == :(x[39] == x[6] * x[7])
		@test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :BILINEAR
		@test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:constr_id] == Set(Any[18])
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[47])]][:y_idx] == 49
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[47])]][:id] == 34
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[47])]][:lifted_constr_ref] == :(x[49] == x[30] * x[47])
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[47])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[47])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[30]), :(x[47])]][:constr_id] == Set(Any[20])
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:y_idx] == 20
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:id] == 5
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:lifted_constr_ref] == :(x[20] == x[5] * x[10])
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:nonlinear_type] == :BINLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:constr_id] == Set(Any[0])
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[6])]][:y_idx] == 21
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[6])]][:id] == 6
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[6])]][:lifted_constr_ref] == :(x[21] == x[11] * x[6])
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[6])]][:nonlinear_type] == :INTLIN
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[6])]][:y_type] == :Cont
		@test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[6])]][:constr_id] == Set(Any[4, 13])
		@test m.internalModel.bounding_constr_mip[1][:rhs] == 3.0
		@test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]
		@test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0, 1.0]
		@test m.internalModel.bounding_constr_mip[1][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[1][:cnt] == 5
		@test m.internalModel.bounding_constr_mip[2][:rhs] == 10.0
		@test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[11]), :(x[12]), :(x[13]), :(x[14]), :(x[15])]
		@test m.internalModel.bounding_constr_mip[2][:coefs] == Any[1.0, 1.0, 1.0, 1.0, 1.0]
		@test m.internalModel.bounding_constr_mip[2][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[2][:cnt] == 5
		@test m.internalModel.bounding_constr_mip[3][:rhs] == 8.0
		@test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[11]), :(x[12]), :(x[13]), :(x[14]), :(x[15])]
		@test m.internalModel.bounding_constr_mip[3][:coefs] == Any[1.0, 1.0, 1.0, 1.0, 1.0]
		@test m.internalModel.bounding_constr_mip[3][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[3][:cnt] == 5
		@test m.internalModel.bounding_constr_mip[4][:rhs] == 17.1
		@test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[21])]
		@test m.internalModel.bounding_constr_mip[4][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[4][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[4][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[5][:rhs] == 17.1
		@test m.internalModel.bounding_constr_mip[5][:vars] == Any[:(x[22])]
		@test m.internalModel.bounding_constr_mip[5][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[5][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[5][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[6][:rhs] == 17.1
		@test m.internalModel.bounding_constr_mip[6][:vars] == Any[:(x[23])]
		@test m.internalModel.bounding_constr_mip[6][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[6][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[6][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[7][:rhs] == 17.1
		@test m.internalModel.bounding_constr_mip[7][:vars] == Any[:(x[24])]
		@test m.internalModel.bounding_constr_mip[7][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[7][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[7][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[8][:rhs] == 17.1
		@test m.internalModel.bounding_constr_mip[8][:vars] == Any[:(x[25])]
		@test m.internalModel.bounding_constr_mip[8][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[8][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[8][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[9][:rhs] == 18.6
		@test m.internalModel.bounding_constr_mip[9][:vars] == Any[:(x[27])]
		@test m.internalModel.bounding_constr_mip[9][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[9][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[9][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[10][:rhs] == 18.6
		@test m.internalModel.bounding_constr_mip[10][:vars] == Any[:(x[29])]
		@test m.internalModel.bounding_constr_mip[10][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[10][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[10][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[11][:rhs] == 18.6
		@test m.internalModel.bounding_constr_mip[11][:vars] == Any[:(x[31])]
		@test m.internalModel.bounding_constr_mip[11][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[11][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[11][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[12][:rhs] == 18.6
		@test m.internalModel.bounding_constr_mip[12][:vars] == Any[:(x[33])]
		@test m.internalModel.bounding_constr_mip[12][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[12][:sense] == :(>=)
		@test m.internalModel.bounding_constr_mip[12][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[13][:rhs] == 28.1
		@test m.internalModel.bounding_constr_mip[13][:vars] == Any[:(x[34])]
		@test m.internalModel.bounding_constr_mip[13][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[13][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[13][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[14][:rhs] == 28.1
		@test m.internalModel.bounding_constr_mip[14][:vars] == Any[:(x[35])]
		@test m.internalModel.bounding_constr_mip[14][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[14][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[14][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[15][:rhs] == 28.1
		@test m.internalModel.bounding_constr_mip[15][:vars] == Any[:(x[36])]
		@test m.internalModel.bounding_constr_mip[15][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[15][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[15][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[16][:rhs] == 28.1
		@test m.internalModel.bounding_constr_mip[16][:vars] == Any[:(x[37])]
		@test m.internalModel.bounding_constr_mip[16][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[16][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[16][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[17][:rhs] == 28.1
		@test m.internalModel.bounding_constr_mip[17][:vars] == Any[:(x[38])]
		@test m.internalModel.bounding_constr_mip[17][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[17][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[17][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[18][:rhs] == 33.2
		@test m.internalModel.bounding_constr_mip[18][:vars] == Any[:(x[42])]
		@test m.internalModel.bounding_constr_mip[18][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[18][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[18][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[19][:rhs] == 33.2
		@test m.internalModel.bounding_constr_mip[19][:vars] == Any[:(x[46])]
		@test m.internalModel.bounding_constr_mip[19][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[19][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[19][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[20][:rhs] == 33.2
		@test m.internalModel.bounding_constr_mip[20][:vars] == Any[:(x[50])]
		@test m.internalModel.bounding_constr_mip[20][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[20][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[20][:cnt] == 1
		@test m.internalModel.bounding_constr_mip[21][:rhs] == 33.2
		@test m.internalModel.bounding_constr_mip[21][:vars] == Any[:(x[54])]
		@test m.internalModel.bounding_constr_mip[21][:coefs] == Any[1.0]
		@test m.internalModel.bounding_constr_mip[21][:sense] == :(<=)
		@test m.internalModel.bounding_constr_mip[21][:cnt] == 1
    end
end

@testset "Expression Parsing || sin/cos" begin
    @testset "Expression Parsing || sin/cos || specialopts " begin

        test_solver=AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)

        m = specialopts(solver=test_solver)
        JuMP.build(m)

        nlk = Vector{Any}(undef, 10)
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
            @test haskey(m.internalModel.nonconvex_terms, k)
            isa(k, Dict) && @test k[:operator] == m.internalModel.nonconvex_terms[k][:nonlinear_type]
        end

        lk = Vector{Any}(undef, 4)
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

    @testset "Expression Parsing || sin/cos || sincos_p1" begin

        test_solver=AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = sincos_p1(solver=test_solver)
        JuMP.build(m)

        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 12
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 2
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[12] == sin(x[1]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 10, 1])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 18
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 8
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[18] == sin(x[4]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 4, 13])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 13
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 3
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[13] == cos(x[2]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 24
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 14
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[24] == sin(x[7]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 7, 16])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 23
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 13
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[23] == cos(x[7]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[16])]][:y_idx] == 33
        @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[16])]][:id] == 23
        @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[16])]][:lifted_constr_ref] == :(x[33] == x[3] * x[4] * x[16])
        @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[16])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[16])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[16])]][:constr_id] == Set(Any[3, 12])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 14
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 4
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[14] == sin(x[2]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 2, 11])
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[18])]][:y_idx] == 34
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[18])]][:id] == 24
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[18])]][:lifted_constr_ref] == :(x[34] == x[4] * x[5] * x[18])
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[18])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[18])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[18])]][:constr_id] == Set(Any[4, 13])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 19
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 9
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[19] == cos(x[5]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 11
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 1
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[11] == cos(x[1]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[1]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 29
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 19
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[29] == cos(x[10]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:y_idx] == 39
        @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:id] == 29
        @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] == :(x[39] == x[9] * x[10] * x[28])
        @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:constr_id] == Set(Any[9, 18])
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[14])]][:y_idx] == 32
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[14])]][:id] == 22
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[14])]][:lifted_constr_ref] == :(x[32] == x[2] * x[3] * x[14])
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[14])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[14])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[14])]][:constr_id] == Set(Any[2, 11])
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[12])]][:y_idx] == 31
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[12])]][:id] == 21
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[12])]][:lifted_constr_ref] == :(x[31] == x[1] * x[2] * x[12])
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[12])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[12])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[12])]][:constr_id] == Set(Any[10, 1])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 30
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 20
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[30] == sin(x[10]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[10]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 26
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 16
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[26] == sin(x[8]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 17, 8])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 22
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 12
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[22] == sin(x[6]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 15, 6])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 28
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 18
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[28] == sin(x[9]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 9, 18])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 20
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 10
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[20] == sin(x[5]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 14, 5])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 17
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 7
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[17] == cos(x[4]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[4]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 16
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 6
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[16] == sin(x[3]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[0, 3, 12])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[22])]][:y_idx] == 36
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[22])]][:id] == 26
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[22])]][:lifted_constr_ref] == :(x[36] == x[6] * x[7] * x[22])
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[22])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[22])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[22])]][:constr_id] == Set(Any[15, 6])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[6]), :(x[20])]][:y_idx] == 35
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[6]), :(x[20])]][:id] == 25
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[6]), :(x[20])]][:lifted_constr_ref] == :(x[35] == x[5] * x[6] * x[20])
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[6]), :(x[20])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[6]), :(x[20])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[6]), :(x[20])]][:constr_id] == Set(Any[14, 5])
        @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9]), :(x[26])]][:y_idx] == 38
        @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9]), :(x[26])]][:id] == 28
        @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9]), :(x[26])]][:lifted_constr_ref] == :(x[38] == x[8] * x[9] * x[26])
        @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9]), :(x[26])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9]), :(x[26])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9]), :(x[26])]][:constr_id] == Set(Any[17, 8])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 25
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 15
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[25] == cos(x[8]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[8]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 15
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 5
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[15] == cos(x[3]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 27
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 17
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[27] == cos(x[9]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 21
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 11
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[21] == cos(x[6]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[6]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[0])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8]), :(x[24])]][:y_idx] == 37
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8]), :(x[24])]][:id] == 27
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8]), :(x[24])]][:lifted_constr_ref] == :(x[37] == x[7] * x[8] * x[24])
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8]), :(x[24])]][:nonlinear_type] == :MULTILINEAR
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8]), :(x[24])]][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8]), :(x[24])]][:constr_id] == Set(Any[7, 16])
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[31])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[2][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[32])]
        @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[3][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[33])]
        @test m.internalModel.bounding_constr_mip[3][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[3][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[3][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[4][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[34])]
        @test m.internalModel.bounding_constr_mip[4][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[4][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[4][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[5][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[5][:vars] == Any[:(x[35])]
        @test m.internalModel.bounding_constr_mip[5][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[5][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[5][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[6][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[6][:vars] == Any[:(x[36])]
        @test m.internalModel.bounding_constr_mip[6][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[6][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[6][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[7][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[7][:vars] == Any[:(x[37])]
        @test m.internalModel.bounding_constr_mip[7][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[7][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[7][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[8][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[8][:vars] == Any[:(x[38])]
        @test m.internalModel.bounding_constr_mip[8][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[8][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[8][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[9][:rhs] == 0.32
        @test m.internalModel.bounding_constr_mip[9][:vars] == Any[:(x[39])]
        @test m.internalModel.bounding_constr_mip[9][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[9][:sense] == :(>=)
        @test m.internalModel.bounding_constr_mip[9][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[10][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[10][:vars] == Any[:(x[31])]
        @test m.internalModel.bounding_constr_mip[10][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[10][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[10][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[11][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[11][:vars] == Any[:(x[32])]
        @test m.internalModel.bounding_constr_mip[11][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[11][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[11][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[12][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[12][:vars] == Any[:(x[33])]
        @test m.internalModel.bounding_constr_mip[12][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[12][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[12][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[13][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[13][:vars] == Any[:(x[34])]
        @test m.internalModel.bounding_constr_mip[13][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[13][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[13][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[14][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[14][:vars] == Any[:(x[35])]
        @test m.internalModel.bounding_constr_mip[14][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[14][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[14][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[15][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[15][:vars] == Any[:(x[36])]
        @test m.internalModel.bounding_constr_mip[15][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[15][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[15][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[16][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[16][:vars] == Any[:(x[37])]
        @test m.internalModel.bounding_constr_mip[16][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[16][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[16][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[17][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[17][:vars] == Any[:(x[38])]
        @test m.internalModel.bounding_constr_mip[17][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[17][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[17][:cnt] == 1
        @test m.internalModel.bounding_constr_mip[18][:rhs] == 0.42
        @test m.internalModel.bounding_constr_mip[18][:vars] == Any[:(x[39])]
        @test m.internalModel.bounding_constr_mip[18][:coefs] == Any[1.0]
        @test m.internalModel.bounding_constr_mip[18][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[18][:cnt] == 1
    end

    @testset "Expression Parsing || sin/cos || trig" begin
        test_solver=AlpineSolver(minlp_solver=pavito_solver,nlp_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = trig(solver=test_solver)
        JuMP.build(m)
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 8
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 3
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[8] == sin(x[7]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[7]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[1])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 10
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 4
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[10] == cos(x[9]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[9]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[1])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 4
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 1
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[4] == sin(x[3]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[3]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[1])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_idx] == 11
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:id] == 5
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:lifted_constr_ref] == :(x[11] == sin(x[2]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:nonlinear_type] == :sin
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[2]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :sin))][:constr_id] == Set(Any[2])
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_idx] == 6
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:id] == 2
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:lifted_constr_ref] == :(x[6] == cos(x[5]))
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:nonlinear_type] == :cos
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:y_type] == :Cont
        @test m.internalModel.nonconvex_terms[Dict{Symbol,Any}(Pair{Symbol,Any}(:vars, Any[5]),Pair{Symbol,Any}(:scalar, 1.0),Pair{Symbol,Any}(:operator, :cos))][:constr_id] == Set(Any[1])
        lk1 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(13.0, 2)])))
        @test m.internalModel.linear_terms[lk1][:y_idx] == 5
        @test m.internalModel.linear_terms[lk1][:id] == 2
        @test m.internalModel.linear_terms[lk1][:y_type] == :(Cont)
        lk2 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(11.0, 2)])))
        @test m.internalModel.linear_terms[lk2][:y_idx] == 3
        @test m.internalModel.linear_terms[lk2][:id] == 1
        @test m.internalModel.linear_terms[lk2][:y_type] == :(Cont)
        lk3 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(19.0, 2)])))
        @test m.internalModel.linear_terms[lk3][:y_idx] == 9
        @test m.internalModel.linear_terms[lk3][:id] == 4
        @test m.internalModel.linear_terms[lk3][:y_type] == :(Cont)
        lk4 = Dict{Symbol,Any}(Pair{Symbol,Any}(:sign, :+),Pair{Symbol,Any}(:scalar, 0.0),Pair{Symbol,Any}(:coef_var, Set(Any[(17.0, 2)])))
        @test m.internalModel.linear_terms[lk4][:y_idx] == 7
        @test m.internalModel.linear_terms[lk4][:id] == 3
        @test m.internalModel.linear_terms[lk4][:y_type] == :(Cont)
        @test m.internalModel.bounding_constr_mip[1][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[4]), :(x[6]), :(x[8]), :(x[10]), :(x[1])]
        @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[-1.0, -1.0, 1.0, 1.0, 1.0]
        @test m.internalModel.bounding_constr_mip[1][:sense] == :(==)
        @test m.internalModel.bounding_constr_mip[1][:cnt] == 5
        @test m.internalModel.bounding_constr_mip[2][:rhs] == 0.0
        @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[11]), :(x[2])]
        @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[5.0, -1.0]
        @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
        @test m.internalModel.bounding_constr_mip[2][:cnt] == 2
    end
end
