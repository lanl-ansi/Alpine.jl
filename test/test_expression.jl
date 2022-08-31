@testset "Expression Parsing || bilinear || Affine || exprs.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = exprstest(solver = test_solver)

    alpine = _build(m)

    ex = alpine.bounding_constr_expr_mip[1]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[1][:coefs]
    @test affdict[:vars] == [:(x[1])]
    @test affdict[:vars] == alpine.bounding_constr_mip[1][:vars]
    @test isapprox(affdict[:rhs], 109.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[1][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[1][:sense]

    ex = alpine.bounding_constr_expr_mip[2]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [3.0, 3.0, 3.0, 3.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[2][:coefs]
    @test affdict[:vars] == [:(x[8]), :(x[9]), :(x[10]), :(x[11])]
    @test affdict[:vars] == alpine.bounding_constr_mip[2][:vars]
    @test isapprox(affdict[:rhs], 111.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[2][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == alpine.bounding_constr_mip[2][:sense]

    ex = alpine.bounding_constr_expr_mip[3]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-1.0, 20.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[3][:coefs]
    @test affdict[:vars] == [:(x[12]), :(x[13])]
    @test affdict[:vars] == alpine.bounding_constr_mip[3][:vars]
    @test isapprox(affdict[:rhs], 222.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[3][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == alpine.bounding_constr_mip[3][:sense]

    # 1.0 * x[12] - 115.0 >= 0.0
    ex = alpine.bounding_constr_expr_mip[4]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[4][:coefs]
    @test affdict[:vars] == [:(x[12])]
    @test affdict[:vars] == alpine.bounding_constr_mip[4][:vars]
    @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[4][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[4][:sense]

    # 1.0 * x[12] - 115.0 <= 0.0
    ex = alpine.bounding_constr_expr_mip[5]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[5][:coefs]
    @test affdict[:vars] == [:(x[12])]
    @test affdict[:vars] == alpine.bounding_constr_mip[5][:vars]
    @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[5][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == alpine.bounding_constr_mip[5][:sense]

    # -1.0 * x[12] - 115.0 >= 0.0
    ex = alpine.bounding_constr_expr_mip[6]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[6][:coefs]
    @test affdict[:vars] == [:(x[12])]
    @test affdict[:vars] == alpine.bounding_constr_mip[6][:vars]
    @test isapprox(affdict[:rhs], 115.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[6][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[6][:sense]

    # (x[1] + 1.0 * x[14]) - 555.0 >= 0.0
    ex = alpine.bounding_constr_expr_mip[7]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0, 1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[7][:coefs]
    @test affdict[:vars] == [:(x[1]), :(x[14])]
    @test affdict[:vars] == alpine.bounding_constr_mip[7][:vars]
    @test isapprox(affdict[:rhs], 555.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[7][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == alpine.bounding_constr_mip[7][:sense]

    # ((x[8] - 7.0 * x[9]) + x[10] + x[4]) - 6666.0 <= 0.0
    ex = alpine.bounding_constr_expr_mip[8]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0, -7.0, 1.0, 1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[8][:coefs]
    @test affdict[:vars] == [:(x[8]), :(x[9]), :(x[10]), :(x[4])]
    @test affdict[:vars] == alpine.bounding_constr_mip[8][:vars]
    @test isapprox(affdict[:rhs], 6666.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[8][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[8][:sense]

    # ((13.0 * x[1] - x[2]) + 30.0 * x[3] + x[4]) - 77.0 >= 0.0
    ex = alpine.bounding_constr_expr_mip[9]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [13.0, -1.0, 30.0, 1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[9][:coefs]
    @test affdict[:vars] == [:(x[1]), :(x[2]), :(x[3]), :(x[4])]
    @test affdict[:vars] == alpine.bounding_constr_mip[9][:vars]
    @test isapprox(affdict[:rhs], 77.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[9][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == alpine.bounding_constr_mip[9][:sense]
end

@testset "Expression Parsing || bilinear || Affine || nlp1.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )
    m = nlp1(solver = test_solver)

    alpine = _build(m)

    ex = alpine.bounding_constr_expr_mip[1]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[1][:coefs]
    @test affdict[:vars] == [:(x[5])]
    @test affdict[:vars] == alpine.bounding_constr_mip[1][:vars]
    @test isapprox(affdict[:rhs], 8.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[1][:rhs]
    @test affdict[:sense] == :(>=)
    @test affdict[:sense] == alpine.bounding_constr_mip[1][:sense]
end

@testset "Expression Parsing || bilinear || Affine || nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = nlp3(solver = test_solver)

    alpine = _build(m)

    ex = alpine.bounding_constr_expr_mip[1]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [0.0025, 0.0025]
    @test affdict[:coefs] == alpine.bounding_constr_mip[1][:coefs]
    @test affdict[:vars] == [:(x[4]), :(x[6])]
    @test affdict[:vars] == alpine.bounding_constr_mip[1][:vars]
    @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[1][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[1][:sense]

    ex = alpine.bounding_constr_expr_mip[2]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-0.0025, 0.0025, 0.0025]
    @test affdict[:coefs] == alpine.bounding_constr_mip[2][:coefs]
    @test affdict[:vars] == [:(x[4]), :(x[5]), :(x[7])]
    @test affdict[:vars] == alpine.bounding_constr_mip[2][:vars]
    @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[2][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[2][:sense]

    ex = alpine.bounding_constr_expr_mip[3]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [-0.01, 0.01]
    @test affdict[:coefs] == alpine.bounding_constr_mip[3][:coefs]
    @test affdict[:vars] == [:(x[5]), :(x[8])]
    @test affdict[:vars] == alpine.bounding_constr_mip[3][:vars]
    @test isapprox(affdict[:rhs], 1.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[3][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[3][:sense]

    ex = alpine.bounding_constr_expr_mip[4]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test (affdict[:coefs] .== [100.0, -1.0, 833.33252]) == [true, true, true]
    @test affdict[:coefs] == alpine.bounding_constr_mip[4][:coefs]
    @test affdict[:vars] == [:(x[1]), :(x[9]), :(x[4])]
    @test affdict[:vars] == alpine.bounding_constr_mip[4][:vars]
    @test isapprox(affdict[:rhs], 83333.333; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[4][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[4][:sense]

    ex = alpine.bounding_constr_expr_mip[5]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0, -1.0, -1250.0, 1250.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[5][:coefs]
    @test affdict[:vars] == [:(x[10]), :(x[11]), :(x[4]), :(x[5])]
    @test affdict[:vars] == alpine.bounding_constr_mip[5][:vars]
    @test isapprox(affdict[:rhs], 0.0; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[5][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[5][:sense]

    ex = alpine.bounding_constr_expr_mip[6]
    affdict = Alpine.expr_linear_to_affine(ex)
    @test affdict[:coefs] == [1.0, -1.0, -2500.0]
    @test affdict[:coefs] == alpine.bounding_constr_mip[6][:coefs]
    @test affdict[:vars] == [:(x[12]), :(x[13]), :(x[5])]
    @test affdict[:vars] == alpine.bounding_constr_mip[6][:vars]
    @test isapprox(affdict[:rhs], -1.25e6; atol = 1e-3)
    @test affdict[:rhs] == alpine.bounding_constr_mip[6][:rhs]
    @test affdict[:sense] == :(<=)
    @test affdict[:sense] == alpine.bounding_constr_mip[6][:sense]
end

@testset "Expression Parsing || bilinear || Simple || bi1.jl " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = operator_c(solver = test_solver)

    alpine = _build(m) # Setup internal model

    @test length(keys(alpine.nonconvex_terms)) == 8
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 1), Expr(:ref, :x, 1)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 2)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 3)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 4), Expr(:ref, :x, 4)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 3)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 4)])

    # TODO setup detailed check on this problem
end
@testset "Expression Parsing || bilinear || Complex || blend029.jl " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => JUNIPER,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = blend029(solver = test_solver)

    alpine = _build(m) # Setup internal model

    @test length(keys(alpine.nonconvex_terms)) == 28
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 37), Expr(:ref, :x, 55)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 38), Expr(:ref, :x, 56)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 37), Expr(:ref, :x, 26)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 43), Expr(:ref, :x, 26)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 47), Expr(:ref, :x, 59)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 48), Expr(:ref, :x, 60)])
    @test haskey(alpine.nonconvex_terms, [Expr(:ref, :x, 47), Expr(:ref, :x, 36)])

    @test alpine.bounding_constr_mip[1][:rhs] == 1.0
    @test alpine.bounding_constr_mip[1][:vars] ==
          Any[:(x[1]), :(x[4]), :(x[7]), :(x[10]), :(x[49])]
    @test alpine.bounding_constr_mip[1][:sense] == :(==)
    @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0, 1.0]
    @test alpine.bounding_constr_mip[1][:cnt] == 5

    @test alpine.bounding_constr_mip[4][:rhs] == 0.1
    @test alpine.bounding_constr_mip[4][:vars] ==
          Any[:(x[4]), :(x[16]), :(x[25]), :(x[34]), :(x[58])]
    @test alpine.bounding_constr_mip[4][:sense] == :(==)
    @test alpine.bounding_constr_mip[4][:coefs] == Any[-1.0, -1.0, -1.0, 1.0, 1.0]
    @test alpine.bounding_constr_mip[4][:cnt] == 5

    @test alpine.bounding_constr_mip[17][:rhs] == -0.14
    @test alpine.bounding_constr_mip[17][:vars] ==
          Any[:(x[11]), :(x[23]), :(x[32]), :(x[64]), :(x[65])]
    @test alpine.bounding_constr_mip[17][:sense] == :(==)
    @test alpine.bounding_constr_mip[17][:coefs] == Any[-1.0, -1.0, -1.0, -1.0, 1.0]
    @test alpine.bounding_constr_mip[17][:cnt] == 5

    @test alpine.bounding_constr_mip[31][:rhs] == 0.0
    @test alpine.bounding_constr_mip[31][:vars] == Any[:(x[13])]
    @test alpine.bounding_constr_mip[31][:sense] == :(>=)
    @test alpine.bounding_constr_mip[31][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[31][:cnt] == 1

    @test alpine.bounding_constr_mip[145][:rhs] == 1.0
    @test alpine.bounding_constr_mip[145][:vars] == Any[:(x[73])]
    @test alpine.bounding_constr_mip[145][:sense] == :(<=)
    @test alpine.bounding_constr_mip[145][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[145][:cnt] == 1

    @test alpine.bounding_constr_mip[67][:rhs] == -1.0
    @test alpine.bounding_constr_mip[67][:vars] == Any[:(x[73])]
    @test alpine.bounding_constr_mip[67][:sense] == :(>=)
    @test alpine.bounding_constr_mip[67][:coefs] == Any[-1.0]
    @test alpine.bounding_constr_mip[67][:cnt] == 1

    @test alpine.bounding_constr_mip[187][:rhs] == 1.0
    @test alpine.bounding_constr_mip[187][:vars] == Any[:(x[79]), :(x[94])]
    @test alpine.bounding_constr_mip[187][:sense] == :(<=)
    @test alpine.bounding_constr_mip[187][:coefs] == Any[1.0, 1.0]
    @test alpine.bounding_constr_mip[187][:cnt] == 2

    @test alpine.bounding_constr_mip[202][:rhs] == 0.04
    @test alpine.bounding_constr_mip[202][:vars] ==
          Any[:(x[103]), :(x[1]), :(x[13]), :(x[25]), :(x[28]), :(x[31])]
    @test alpine.bounding_constr_mip[202][:sense] == :(==)
    @test alpine.bounding_constr_mip[202][:coefs] == Any[1.0, -0.6, -0.2, 0.2, 0.2, 0.2]
    @test alpine.bounding_constr_mip[202][:cnt] == 6

    @test alpine.nonconvex_terms[[:(x[37]), :(x[55])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[37]), :(x[55])]][:lifted_var_ref] == :(x[103])
    @test alpine.nonconvex_terms[[:(x[37]), :(x[55])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[37]), :(x[55])]][:nonlinear_type] == :BILINEAR

    @test alpine.bounding_constr_mip[206][:rhs] == 0.0
    @test alpine.bounding_constr_mip[206][:vars] ==
          Any[:(x[107]), :(x[103]), :(x[108]), :(x[109]), :(x[110]), :(x[2]), :(x[14])]
    @test alpine.bounding_constr_mip[206][:sense] == :(==)
    @test alpine.bounding_constr_mip[206][:coefs] ==
          Any[1.0, -1.0, 1.0, 1.0, 1.0, -0.6, -0.2]
    @test alpine.bounding_constr_mip[206][:cnt] == 7

    @test alpine.nonconvex_terms[[:(x[38]), :(x[56])]][:id] == 5
    @test alpine.nonconvex_terms[[:(x[38]), :(x[56])]][:lifted_var_ref] == :(x[107])
    @test alpine.nonconvex_terms[[:(x[38]), :(x[56])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[38]), :(x[56])]][:nonlinear_type] == :BILINEAR

    @test alpine.nonconvex_terms[[:(x[37]), :(x[26])]][:id] == 6
    @test alpine.nonconvex_terms[[:(x[37]), :(x[26])]][:lifted_var_ref] == :(x[108])
    @test alpine.nonconvex_terms[[:(x[37]), :(x[26])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[37]), :(x[26])]][:nonlinear_type] == :BILINEAR

    @test alpine.nonconvex_terms[[:(x[37]), :(x[32])]][:id] == 8
    @test alpine.nonconvex_terms[[:(x[37]), :(x[32])]][:lifted_var_ref] == :(x[110])
    @test alpine.nonconvex_terms[[:(x[37]), :(x[32])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[37]), :(x[32])]][:nonlinear_type] == :BILINEAR

    @test alpine.bounding_constr_mip[213][:rhs] == 0.0
    @test alpine.bounding_constr_mip[213][:vars] ==
          Any[:(x[129]), :(x[127]), :(x[124]), :(x[130]), :(x[6]), :(x[18])]
    @test alpine.bounding_constr_mip[213][:sense] == :(==)
    @test alpine.bounding_constr_mip[213][:coefs] == Any[1.0, -1.0, -1.0, 1.0, -0.4, -0.4]
    @test alpine.bounding_constr_mip[213][:cnt] == 6

    @test alpine.nonconvex_terms[[:(x[48]), :(x[60])]][:id] == 27
    @test alpine.nonconvex_terms[[:(x[48]), :(x[60])]][:lifted_var_ref] == :(x[129])
    @test alpine.nonconvex_terms[[:(x[48]), :(x[60])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[48]), :(x[60])]][:nonlinear_type] == :BILINEAR

    @test alpine.nonconvex_terms[[:(x[47]), :(x[59])]][:id] == 25
    @test alpine.nonconvex_terms[[:(x[47]), :(x[59])]][:lifted_var_ref] == :(x[127])
    @test alpine.nonconvex_terms[[:(x[47]), :(x[59])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[47]), :(x[59])]][:nonlinear_type] == :BILINEAR

    @test alpine.nonconvex_terms[[:(x[47]), :(x[36])]][:id] == 28
    @test alpine.nonconvex_terms[[:(x[47]), :(x[36])]][:lifted_var_ref] == :(x[130])
    @test alpine.nonconvex_terms[[:(x[47]), :(x[36])]][:convexified] == false
    @test alpine.nonconvex_terms[[:(x[47]), :(x[36])]][:nonlinear_type] == :BILINEAR
end

@testset "Expression Parsing || multilinear || Simple || multi.jl " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = multi3(solver = test_solver, exprmode = 1)

    alpine = _build(m) # Setup internal model

    @test length(keys(alpine.nonconvex_terms)) == 1
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:lifted_var_ref] == :(x[4])
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] ==
          :MULTILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[4])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    @test alpine.bounding_constr_mip[1][:rhs] == 3.0
    @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3])]
    @test alpine.bounding_constr_mip[1][:sense] == :(<=)
    @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0]
    @test alpine.bounding_constr_mip[1][:cnt] == 3

    m = multi3(solver = test_solver, exprmode = 2)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] ==
          :(x[4])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 3)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[5])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi3(solver = test_solver, exprmode = 3)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] ==
          :(x[4])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[5])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 1)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 1

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[5])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] ==
          :MULTILINEAR

    @test alpine.bounding_constr_mip[1][:rhs] == 4.0
    @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2]), :(x[3]), :(x[4])]
    @test alpine.bounding_constr_mip[1][:sense] == :(<=)
    @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0, 1.0, 1.0]
    @test alpine.bounding_constr_mip[1][:cnt] == 4

    m = multi4(solver = test_solver, exprmode = 2)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 3

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:id] == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:lifted_var_ref] ==
          :(x[7])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 6)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 3)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[:(x[5]), :(x[3]), :(x[4])]][:id] == 2
    @test alpine.nonconvex_terms[[:(x[5]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[6])
    @test alpine.nonconvex_terms[[:(x[5]), :(x[3]), :(x[4])]][:nonlinear_type] ==
          :MULTILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 4)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[5])]][:id] == 2
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[5])]][:lifted_var_ref] == :(x[6])
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[5])]][:nonlinear_type] ==
          :MULTILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 5)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 3

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 5)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[7])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 6)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2

    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:lifted_var_ref] == :(x[5])
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 5), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 7)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 3

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 3), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 5)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] ==
          :(x[7])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 8)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[:(x[1]), :(x[5]), :(x[4])]][:id] == 2
    @test alpine.nonconvex_terms[[:(x[1]), :(x[5]), :(x[4])]][:lifted_var_ref] == :(x[6])
    @test alpine.nonconvex_terms[[:(x[1]), :(x[5]), :(x[4])]][:nonlinear_type] ==
          :MULTILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 9)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 2

    @test alpine.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4])]][:lifted_var_ref] == :(x[5])
    @test alpine.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[6])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 10)

    alpine = _build(m)
    @test length(keys(alpine.nonconvex_terms)) == 3

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 4), Expr(:ref, :x, 5)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:id] == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:lifted_var_ref] ==
          :(x[7])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 6)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing

    m = multi4(solver = test_solver, exprmode = 11)

    alpine = _build(m)

    @test length(keys(alpine.nonconvex_terms)) == 3

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 3)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:lifted_var_ref] ==
          :(x[6])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 5)]][:nonlinear_type] ==
          :BILINEAR
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:id] == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:lifted_var_ref] ==
          :(x[7])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 6), Expr(:ref, :x, 4)]][:nonlinear_type] ==
          :BILINEAR

    @test alpine.bounding_obj_mip[:rhs] == 0
    @test alpine.bounding_obj_mip[:vars] == Expr[:(x[7])]
    @test alpine.bounding_obj_mip[:coefs] == [1.0]
    @test alpine.bounding_obj_mip[:cnt] == 1
    @test alpine.bounding_obj_mip[:sense] === nothing
end

@testset "Expression Parsing || bilinear || Complex-div || div.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = div(solver = test_solver)

    alpine = _build(m) # Setup internal model

    @test length(keys(alpine.nonconvex_terms)) == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:id] == 1
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:lifted_var_ref] ==
          :(x[3])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 1)]][:nonlinear_type] ==
          :MONOMIAL

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:id] == 2
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:lifted_var_ref] ==
          :(x[4])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 2), Expr(:ref, :x, 2)]][:nonlinear_type] ==
          :MONOMIAL

    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:id] == 3
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:lifted_var_ref] ==
          :(x[5])
    @test alpine.nonconvex_terms[[Expr(:ref, :x, 1), Expr(:ref, :x, 2)]][:nonlinear_type] ==
          :BILINEAR

    aff_mip = alpine.bounding_constr_mip

    @test aff_mip[1][:rhs] == 0.0
    @test aff_mip[1][:vars] == Any[:(x[1])]
    @test aff_mip[1][:sense] == :(>=)
    @test round(aff_mip[1][:coefs][1]; digits = 1) == -3.0
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
    @test aff_mip[8][:vars] == Any[:(x[2]), :(x[5])]
    @test aff_mip[8][:sense] == :(>=)
    @test aff_mip[8][:coefs] == Any[5.6, -72000.0]
    @test aff_mip[8][:cnt] == 2

    @test aff_mip[9][:rhs] == 0.0
    @test aff_mip[9][:vars] == Any[:(x[2]), :(x[5])]
    @test aff_mip[9][:sense] == :(>=)
    @test aff_mip[9][:coefs] == Any[5.6, -36000.0]
    @test aff_mip[9][:cnt] == 2

    @test aff_mip[10][:rhs] == 0.0
    @test aff_mip[10][:vars] == Any[:(x[2]), :(x[5]), :(x[2])]
    @test aff_mip[10][:sense] == :(>=)
    @test aff_mip[10][:coefs] == Any[5.6, -300.0, -1.75]
    @test aff_mip[10][:cnt] == 3

    @test aff_mip[11][:rhs] == 0.0
    @test aff_mip[11][:vars] == Any[:(x[2]), :(x[1]), :(x[5])]
    @test aff_mip[11][:sense] == :(>=)
    @test aff_mip[11][:coefs] == Any[5.6, -0.5, 0.5]
    @test aff_mip[11][:cnt] == 3
end

@testset "Expression Parsing || part1 " begin
    m = Model(
        optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        ),
    )
    @variable(m, x[1:4] >= 0)
    @NLconstraint(m, x[1]^2 >= 1)    # Basic monomial x[5]=x[1]^2
    @NLconstraint(m, x[1] * x[2] <= 1)  # x[6] <= 1 : x[6] = x[1]*x[2]
    @NLconstraint(m, x[1]^2 * x[2]^2 <= 1)          # x[5] + x[7] <= 1 : x[7] = x[2]^2
    @NLconstraint(m, x[1] * (x[2] * x[3]) >= 1)         # x[9] >= 1 : x[8] = x[2] * x[3] && x[9] = x[1]*x[8]
    @NLconstraint(m, x[1]^2 * (x[2]^2 * x[3]^2) <= 1) # x[12] <= 1 : x[10] = x[3] ^ 2 && x[11] = x[7] * x[10] && x[12] = x[11]*[5]

    alpine = _build(m)

    @test alpine.bounding_constr_expr_mip[1] == :(x[5] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[2] == :(x[6] - 1.0 <= 0.0)
    @test alpine.bounding_constr_expr_mip[3] == :(x[8] - 1.0 <= 0.0)
    @test alpine.bounding_constr_expr_mip[4] == :(x[10] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[5] == :(x[13] - 1.0 <= 0.0)
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[1])])
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2])])
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[2])])
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[3])])
    @test haskey(alpine.nonconvex_terms, [:(x[5]), :(x[7])])
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[9])])
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[3])])
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2])])
    @test haskey(alpine.nonconvex_terms, [:(x[5]), :(x[12])])

    @test alpine.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2])]][:id] == 2
    @test alpine.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 3
    @test alpine.nonconvex_terms[[:(x[5]), :(x[7])]][:id] == 4
    @test alpine.nonconvex_terms[[:(x[2]), :(x[3])]][:id] == 5
    @test alpine.nonconvex_terms[[:(x[1]), :(x[9])]][:id] == 6
    @test alpine.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 7
    @test alpine.nonconvex_terms[[:(x[7]), :(x[11])]][:id] == 8
    @test alpine.nonconvex_terms[[:(x[5]), :(x[12])]][:id] == 9
end

@testset "Expression Parsing || part2" begin
    m = Model(
        optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        ),
    )

    @variable(m, x[1:4] >= 0)
    @NLconstraint(m, (x[1] * x[2]) * x[3] >= 1)
    @NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

    @NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
    @NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2 <= 1)
    @NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
    @NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

    alpine = _build(m)

    @test alpine.bounding_constr_expr_mip[1] == :(x[6] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[2] == :(x[11] - 1.0 <= 0.0)
    @test alpine.bounding_constr_expr_mip[3] == :(x[13] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[4] == :(x[15] - 1.0 <= 0.0)
    @test alpine.bounding_constr_expr_mip[5] == :(x[17] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[6] == :(x[19] - 1.0 <= 0.0)

    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2])]) #5
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[5])]) #6
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[1])]) #7
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[2])]) #8
    @test haskey(alpine.nonconvex_terms, [:(x[7]), :(x[8])]) #9
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[3])]) #10
    @test haskey(alpine.nonconvex_terms, [:(x[9]), :(x[10])]) #11
    @test haskey(alpine.nonconvex_terms, [:(x[8]), :(x[10])]) #12
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[12])]) #13
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[7])])  #14
    @test haskey(alpine.nonconvex_terms, [:(x[14]), :(x[10])]) #15
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[3])]) #16
    @test haskey(alpine.nonconvex_terms, [:(x[7]), :(x[16])]) #17
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[8])]) #18
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[18])]) #19

    @test alpine.nonconvex_terms[[:(x[1]), :(x[2])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[3]), :(x[5])]][:id] == 2
    @test alpine.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 3
    @test alpine.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 4
    @test alpine.nonconvex_terms[[:(x[7]), :(x[8])]][:id] == 5
    @test alpine.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 6
    @test alpine.nonconvex_terms[[:(x[9]), :(x[10])]][:id] == 7
    @test alpine.nonconvex_terms[[:(x[8]), :(x[10])]][:id] == 8
    @test alpine.nonconvex_terms[[:(x[1]), :(x[12])]][:id] == 9
    @test alpine.nonconvex_terms[[:(x[2]), :(x[7])]][:id] == 10
    @test alpine.nonconvex_terms[[:(x[14]), :(x[10])]][:id] == 11
    @test alpine.nonconvex_terms[[:(x[2]), :(x[3])]][:id] == 12
    @test alpine.nonconvex_terms[[:(x[7]), :(x[16])]][:id] == 13
    @test alpine.nonconvex_terms[[:(x[1]), :(x[8])]][:id] == 14
    @test alpine.nonconvex_terms[[:(x[3]), :(x[18])]][:id] == 15
end

@testset "Expression Parsing || part3" begin
    m = Model(
        optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        ),
    )

    @variable(m, x[1:4] >= 0)
    @NLconstraint(m, ((x[1] * x[2]) * x[3]) * x[4] >= 1)
    @NLconstraint(m, ((x[1]^2 * x[2]) * x[3]) * x[4] <= 1)
    @NLconstraint(m, ((x[1] * x[2]^2) * x[3]) * x[4] >= 1)
    @NLconstraint(m, ((x[1] * x[2]) * x[3]^2) * x[4] <= 1)
    @NLconstraint(m, ((x[1] * x[2]) * x[3]) * x[4]^2 >= 1)
    @NLconstraint(m, ((x[1]^2 * x[2]^2) * x[3]^2) * x[4]^2 <= 1)

    alpine = _build(m)

    @test alpine.bounding_constr_expr_mip[1] == :(x[7] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[2] == :(x[11] - 1.0 <= 0.0)
    @test alpine.bounding_constr_expr_mip[3] == :(x[15] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[4] == :(x[18] - 1.0 <= 0.0)
    @test alpine.bounding_constr_expr_mip[5] == :(x[20] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[6] == :(x[23] - 1.0 <= 0.0)

    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2])]) #5
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[5])]) #6
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[6])]) #7
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[1])]) #8
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[8])]) #9
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[9])]) #10
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[10])]) #11
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[2])]) #12
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[12])])  #13
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[13])]) #14
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[14])]) #15
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[3])]) #16
    @test haskey(alpine.nonconvex_terms, [:(x[5]), :(x[16])]) #17
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[17])]) #18
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[4])]) #19
    @test haskey(alpine.nonconvex_terms, [:(x[6]), :(x[19])]) #20
    @test haskey(alpine.nonconvex_terms, [:(x[8]), :(x[12])]) #21
    @test haskey(alpine.nonconvex_terms, [:(x[21]), :(x[16])]) #22
    @test haskey(alpine.nonconvex_terms, [:(x[22]), :(x[19])]) #23

    @test alpine.nonconvex_terms[[:(x[1]), :(x[2])]][:id] == 1
    @test alpine.nonconvex_terms[[:(x[3]), :(x[5])]][:id] == 2
    @test alpine.nonconvex_terms[[:(x[4]), :(x[6])]][:id] == 3
    @test alpine.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 4
    @test alpine.nonconvex_terms[[:(x[2]), :(x[8])]][:id] == 5
    @test alpine.nonconvex_terms[[:(x[3]), :(x[9])]][:id] == 6
    @test alpine.nonconvex_terms[[:(x[4]), :(x[10])]][:id] == 7
    @test alpine.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 8
    @test alpine.nonconvex_terms[[:(x[1]), :(x[12])]][:id] == 9
    @test alpine.nonconvex_terms[[:(x[3]), :(x[13])]][:id] == 10
    @test alpine.nonconvex_terms[[:(x[4]), :(x[14])]][:id] == 11
    @test alpine.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 12
    @test alpine.nonconvex_terms[[:(x[5]), :(x[16])]][:id] == 13
    @test alpine.nonconvex_terms[[:(x[4]), :(x[17])]][:id] == 14
    @test alpine.nonconvex_terms[[:(x[4]), :(x[4])]][:id] == 15
    @test alpine.nonconvex_terms[[:(x[6]), :(x[19])]][:id] == 16
    @test alpine.nonconvex_terms[[:(x[8]), :(x[12])]][:id] == 17
    @test alpine.nonconvex_terms[[:(x[21]), :(x[16])]][:id] == 18
    @test alpine.nonconvex_terms[[:(x[22]), :(x[19])]][:id] == 19
end

@testset "Expression Parsing || part7" begin
    m = Model(
        optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        ),
    )
    @variable(m, x[1:4] >= 0)

    @NLconstraint(m, x[1] * x[2] * x[3] * x[4] >= 1)
    @NLconstraint(m, x[1]^2 * x[2] * x[3] * x[4] >= 1)
    @NLconstraint(m, x[1] * x[2]^2 * x[3] * x[4] >= 1)
    @NLconstraint(m, x[1] * x[2] * x[3]^2 * x[4]^2 >= 1)
    @NLconstraint(m, x[1] * x[2]^2 * x[3]^2 * x[4] >= 1)
    @NLconstraint(m, x[1]^2 * x[2] * x[3] * x[4]^2 >= 1)
    @NLconstraint(m, x[1]^2 * x[2]^2 * x[3]^2 * x[4]^2 >= 1)

    alpine = _build(m)

    @test alpine.bounding_constr_expr_mip[1] == :(x[5] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[2] == :(x[7] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[3] == :(x[9] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[4] == :(x[12] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[5] == :(x[13] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[6] == :(x[14] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[7] == :(x[15] - 1.0 >= 0.0)

    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2]), :(x[3]), :(x[4])]) #5
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[1])]) #6
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[3]), :(x[4]), :(x[6])]) #7
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[2])]) #8
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[3]), :(x[4]), :(x[8])]) #9
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[3])]) #10
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[4])]) #11
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2]), :(x[10]), :(x[11])]) #12
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[4]), :(x[8]), :(x[10])]) #13
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[3]), :(x[6]), :(x[11])]) #14
    @test haskey(alpine.nonconvex_terms, [:(x[6]), :(x[8]), :(x[10]), :(x[11])]) #15

    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 1 #5
    @test alpine.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 2 #6
    @test alpine.nonconvex_terms[[:(x[2]), :(x[3]), :(x[4]), :(x[6])]][:id] == 3 #7
    @test alpine.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 4 #8
    @test alpine.nonconvex_terms[[:(x[1]), :(x[3]), :(x[4]), :(x[8])]][:id] == 5 #9
    @test alpine.nonconvex_terms[[:(x[3]), :(x[3])]][:id] == 6 #10
    @test alpine.nonconvex_terms[[:(x[4]), :(x[4])]][:id] == 7 #11
    @test alpine.nonconvex_terms[[:(x[1]), :(x[2]), :(x[10]), :(x[11])]][:id] == 8  #12
    @test alpine.nonconvex_terms[[:(x[1]), :(x[4]), :(x[8]), :(x[10])]][:id] == 9 #13
    @test alpine.nonconvex_terms[[:(x[2]), :(x[3]), :(x[6]), :(x[11])]][:id] == 10 #14
    @test alpine.nonconvex_terms[[:(x[6]), :(x[8]), :(x[10]), :(x[11])]][:id] == 11 #15
end

@testset "Expression Parsing || part8" begin
    m = Model(
        optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        ),
    )
    @variable(m, x[1:4] >= 0)

    @NLconstraint(m, (x[1] * x[2] * x[3]) * x[4] >= 1)
    @NLconstraint(m, (x[1]^2 * x[2] * x[3]) * x[4] >= 1)
    @NLconstraint(m, x[1] * (x[2]^2 * x[3]) * x[4] >= 1)
    @NLconstraint(m, x[1] * (x[2] * x[3]^2) * x[4] >= 1)
    @NLconstraint(m, (x[1] * x[2]^2) * x[3] * x[4] >= 1)
    @NLconstraint(m, (x[1] * x[2]) * x[3]^2 * x[4] >= 1)
    @NLconstraint(m, (x[1] * x[2]) * x[3] * x[4]^2 >= 1)
    @NLconstraint(m, (x[1] * x[2]) * x[3]^2 * x[4]^2 >= 1)
    @NLconstraint(m, (x[1]^2 * x[2]^2 * x[3]^2) * x[4]^2 >= 1)

    alpine = _build(m)

    @test alpine.bounding_constr_expr_mip[1] == :(x[6] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[2] == :(x[9] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[3] == :(x[12] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[4] == :(x[15] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[5] == :(x[17] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[6] == :(x[19] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[7] == :(x[21] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[8] == :(x[22] - 1.0 >= 0.0)
    @test alpine.bounding_constr_expr_mip[9] == :(x[24] - 1.0 >= 0.0)

    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2]), :(x[3])]) #5
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[5])]) #6
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[1])]) #7
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[3]), :(x[7])]) #8
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[8])]) #9
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[2])]) #10
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[10])]) #11
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[4]), :(x[11])]) #12
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[3])]) #13
    @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[13])]) #14
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[4]), :(x[14])]) #15
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[10])]) #16
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[4]), :(x[16])]) #17
    @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[2])]) #18
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[18]), :(x[13])]) #19
    @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[4])]) #20
    @test haskey(alpine.nonconvex_terms, [:(x[3]), :(x[18]), :(x[20])]) #21
    @test haskey(alpine.nonconvex_terms, [:(x[18]), :(x[13]), :(x[20])]) #22
    @test haskey(alpine.nonconvex_terms, [:(x[7]), :(x[10]), :(x[13])]) #23
    @test haskey(alpine.nonconvex_terms, [:(x[23]), :(x[20])]) #24
end

@testset "Expression Parsing || Convex" begin
    @testset "Convex Parsing :: PART I" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        )
        m = convex_test(solver = test_solver)

        alpine = _build(m)

        @test alpine.num_constr_convex == 21

        # 0 : OBJ
        @test alpine.obj_structure == :convex
        @test alpine.nonlinear_constrs[0][:expr_orig] == :objective
        @test alpine.nonlinear_constrs[0][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[0][:convexified] == false

        @test alpine.bounding_obj_mip[:sense] === nothing
        @test alpine.bounding_obj_mip[:coefs] == [1.0, 1.0]
        @test alpine.bounding_obj_mip[:vars] == [:(x[1]), :(x[3])]
        @test alpine.bounding_obj_mip[:rhs] == 0.0
        @test alpine.bounding_obj_mip[:powers] == [2, 2]
        @test alpine.bounding_obj_mip[:cnt] == 2

        # 1
        @test alpine.constr_structure[1] == :convex
        @test alpine.nonlinear_constrs[1][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[1][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[1][:convexified] == false
        @test alpine.bounding_constr_mip[1][:sense] == :(<=)
        @test alpine.bounding_constr_mip[1][:coefs] == [3.0, 4.0]
        @test alpine.bounding_constr_mip[1][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[1][:rhs] == 25.0
        @test alpine.bounding_constr_mip[1][:powers] == [2, 2]
        @test alpine.bounding_constr_mip[1][:cnt] == 2

        # 2
        @test alpine.constr_structure[2] == :convex
        @test alpine.nonlinear_constrs[2][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[2][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[2][:convexified] == false
        @test alpine.bounding_constr_mip[2][:sense] == :(<=)
        @test alpine.bounding_constr_mip[2][:coefs] == [3.0, 4.0]
        @test alpine.bounding_constr_mip[2][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[2][:rhs] == 25.0
        @test alpine.bounding_constr_mip[2][:powers] == [2, 2]
        @test alpine.bounding_constr_mip[2][:cnt] == 2

        # 4
        @test alpine.constr_structure[4] == :convex
        @test alpine.nonlinear_constrs[4][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[4][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[4][:convexified] == false
        @test alpine.bounding_constr_mip[4][:sense] == :(<=)
        @test alpine.bounding_constr_mip[4][:coefs] == [3.0, 4.0]
        @test alpine.bounding_constr_mip[4][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[4][:rhs] == 10.0
        @test alpine.bounding_constr_mip[4][:powers] == [2, 2]
        @test alpine.bounding_constr_mip[4][:cnt] == 2

        # 5
        @test alpine.constr_structure[5] == :convex
        @test alpine.nonlinear_constrs[5][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[5][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[5][:convexified] == :false
        @test alpine.bounding_constr_mip[5][:sense] == :(<=)
        @test alpine.bounding_constr_mip[5][:coefs] == [3.0, 4.0, 6.0]
        @test alpine.bounding_constr_mip[5][:vars] == [:(x[1]), :(x[2]), :(x[3])]
        @test alpine.bounding_constr_mip[5][:rhs] == 10.0
        @test alpine.bounding_constr_mip[5][:powers] == [2, 2, 2]
        @test alpine.bounding_constr_mip[5][:cnt] == 3

        # 6
        @test alpine.constr_structure[6] == :convex
        @test alpine.nonlinear_constrs[6][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[6][:convex_type] == :convexC
        @test alpine.nonlinear_constrs[6][:convexified] == :false
        @test alpine.bounding_constr_mip[6][:sense] == :(<=)
        @test alpine.bounding_constr_mip[6][:coefs] == [3.0, 4.0, 5.0]
        @test alpine.bounding_constr_mip[6][:vars] == [:(x[1]), :(x[2]), :(x[5])]
        @test alpine.bounding_constr_mip[6][:rhs] == 100.0
        @test alpine.bounding_constr_mip[6][:powers] == [0.5, 0.5, 0.5]
        @test alpine.bounding_constr_mip[6][:cnt] == 3

        # 7
        @test alpine.constr_structure[7] == :convex
        @test alpine.nonlinear_constrs[7][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[7][:convex_type] == :convexC
        @test alpine.nonlinear_constrs[7][:convexified] == :false
        @test alpine.bounding_constr_mip[7][:sense] == :(>=)
        @test alpine.bounding_constr_mip[7][:coefs] == [-3.0, -4.0]
        @test alpine.bounding_constr_mip[7][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[7][:rhs] == -100.0
        @test alpine.bounding_constr_mip[7][:powers] == [0.5, 0.5]
        @test alpine.bounding_constr_mip[7][:cnt] == 2

        # 8
        @test alpine.constr_structure[8] == :convex
        @test alpine.nonlinear_constrs[8][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[8][:convex_type] == :convexB
        @test alpine.nonlinear_constrs[8][:convexified] == :false
        @test alpine.bounding_constr_mip[8][:sense] == :(<=)
        @test alpine.bounding_constr_mip[8][:coefs] == [3.0, 1.0, 5.0]
        @test alpine.bounding_constr_mip[8][:vars] == [:(x[1]), :(x[2]), :(x[3])]
        @test alpine.bounding_constr_mip[8][:rhs] == 200.0
        @test alpine.bounding_constr_mip[8][:powers] == [3.0, 3.0, 3.0]
        @test alpine.bounding_constr_mip[8][:cnt] == 3

        # 9
        @test alpine.constr_structure[9] == :convex
        @test alpine.nonlinear_constrs[9][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[9][:convex_type] == :convexB
        @test alpine.nonlinear_constrs[9][:convexified] == :false
        @test alpine.bounding_constr_mip[9][:sense] == :(<=)
        @test alpine.bounding_constr_mip[9][:coefs] == [1.0, 1.0, 1.0, 100.0]
        @test alpine.bounding_constr_mip[9][:vars] == [:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        @test alpine.bounding_constr_mip[9][:rhs] == 200.0
        @test alpine.bounding_constr_mip[9][:powers] == [3, 3, 3, 3]
        @test alpine.bounding_constr_mip[9][:cnt] == 4

        # 11
        @test alpine.constr_structure[11] == :convex
        @test alpine.nonlinear_constrs[11][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[11][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[11][:convexified] == :false
        @test alpine.bounding_constr_mip[11][:sense] == :(<=)
        @test alpine.bounding_constr_mip[11][:coefs] == [3.0, 4.0]
        @test alpine.bounding_constr_mip[11][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[11][:rhs] == 25.0
        @test alpine.bounding_constr_mip[11][:powers] == [2, 2]
        @test alpine.bounding_constr_mip[11][:cnt] == 2

        # 14
        @test alpine.constr_structure[14] == :convex
        @test alpine.nonlinear_constrs[14][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[14][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[14][:convexified] == :false
        @test alpine.bounding_constr_mip[14][:sense] == :(<=)
        @test alpine.bounding_constr_mip[14][:coefs] == [3.0, 5.0]
        @test alpine.bounding_constr_mip[14][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[14][:rhs] == 25.0
        @test alpine.bounding_constr_mip[14][:powers] == [2, 2]
        @test alpine.bounding_constr_mip[14][:cnt] == 2

        # 15
        @test alpine.constr_structure[15] == :convex
        @test alpine.nonlinear_constrs[15][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[15][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[15][:convexified] == :false
        @test alpine.bounding_constr_mip[15][:sense] == :(<=)
        @test alpine.bounding_constr_mip[15][:coefs] == [3.0, 5.0, 1.0]
        @test alpine.bounding_constr_mip[15][:vars] == [:(x[1]), :(x[2]), :(x[4])]
        @test alpine.bounding_constr_mip[15][:rhs] == 25.0
        @test alpine.bounding_constr_mip[15][:powers] == [2, 2, 2]
        @test alpine.bounding_constr_mip[15][:cnt] == 3

        # 19
        @test alpine.constr_structure[19] == :convex
        @test alpine.nonlinear_constrs[19][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[19][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[19][:convexified] == :false
        @test alpine.bounding_constr_mip[19][:sense] == :(<=)
        @test alpine.bounding_constr_mip[19][:coefs] == [3.0, 16.0]
        @test alpine.bounding_constr_mip[19][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[19][:rhs] == 40.0
        @test alpine.bounding_constr_mip[19][:powers] == [2, 2]
        @test alpine.bounding_constr_mip[19][:cnt] == 2

        # 22
        @test alpine.constr_structure[22] == :convex
        @test alpine.nonlinear_constrs[22][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[22][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[22][:convexified] == :false
        @test alpine.bounding_constr_mip[22][:sense] == :(<=)
        @test alpine.bounding_constr_mip[22][:coefs] == [3.0, 4.0, 5.0, 6.0]
        @test alpine.bounding_constr_mip[22][:vars] ==
              [:(x[1]), :(x[2]), :(x[3]), :(x[4])]
        @test alpine.bounding_constr_mip[22][:rhs] == 15.0
        @test alpine.bounding_constr_mip[22][:powers] == [2, 2, 2, 2]
        @test alpine.bounding_constr_mip[22][:cnt] == 4

        # 25
        @test alpine.constr_structure[25] == :convex
        @test alpine.nonlinear_constrs[25][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[25][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[25][:convexified] == :false
        @test alpine.bounding_constr_mip[25][:sense] == :(<=)
        @test alpine.bounding_constr_mip[25][:coefs] == [1.0, 1.0, 1.0, 1.0, 1.0]
        @test alpine.bounding_constr_mip[25][:vars] ==
              [:(x[1]), :(x[2]), :(x[3]), :(x[4]), :(x[5])]
        @test alpine.bounding_constr_mip[25][:rhs] == 99999.0
        @test alpine.bounding_constr_mip[25][:powers] == [2, 2, 2, 2, 2]
        @test alpine.bounding_constr_mip[25][:cnt] == 5

        # 26
        @test alpine.constr_structure[26] == :convex
        @test alpine.nonlinear_constrs[26][:expr_orig] == :constraints
        @test alpine.nonlinear_constrs[26][:convex_type] == :convexA
        @test alpine.nonlinear_constrs[26][:convexified] == :false
        @test alpine.bounding_constr_mip[26][:sense] == :(<=)
        @test alpine.bounding_constr_mip[26][:coefs] == [3.0, 4.0]
        @test alpine.bounding_constr_mip[26][:vars] == [:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[26][:rhs] == 200.0
        @test alpine.bounding_constr_mip[26][:powers] == [4, 4]
        @test alpine.bounding_constr_mip[26][:cnt] == 2
    end
end

@testset "Expression Prasing || Linear Lifting" begin
    @testset "Expression Parsing || Linear Lifting || nlp2" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "disc_ratio" => 8,
            "log_level" => 100,
        )

        m = nlp2(solver = test_solver)

        alpine = _build(m)

        @test length(alpine.linear_terms) == 2
        @test length(alpine.nonconvex_terms) == 4

        lk = Vector{Any}(undef, 2)
        lk[1] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, -1.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3)])),
        )
        lk[2] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, -2.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6)])),
        )

        @test length(keys(alpine.linear_terms)) == 2
        yids = [4, 7]
        for i in 1:length(keys(alpine.linear_terms))
            for j in keys(alpine.linear_terms)
                if alpine.linear_terms[j][:id] == i
                    @test j == lk[i]
                    @test j[:sign] == :+
                    @test j[:scalar] == lk[i][:scalar]
                    @test j[:coef_var] == lk[i][:coef_var]
                    @test alpine.linear_terms[j][:y_idx] == yids[i]
                end
            end
        end

        @test haskey(alpine.linear_terms, lk[1])
        @test haskey(alpine.linear_terms, lk[2])

        @test haskey(alpine.nonconvex_terms, [:(x[2]), :(x[2])])
        @test haskey(alpine.nonconvex_terms, [:(x[4]), :(x[4])])
        @test haskey(alpine.nonconvex_terms, [:(x[7]), :(x[7])])
        @test haskey(alpine.nonconvex_terms, [:(x[1]), :(x[1])])
        @test alpine.nonconvex_terms[[:(x[1]), :(x[1])]][:id] == 1
        @test alpine.nonconvex_terms[[:(x[2]), :(x[2])]][:id] == 3
        @test alpine.nonconvex_terms[[:(x[4]), :(x[4])]][:id] == 2
        @test alpine.nonconvex_terms[[:(x[7]), :(x[7])]][:id] == 4
        @test alpine.nonconvex_terms[[:(x[1]), :(x[1])]][:lifted_var_ref].args[2] == 3
        @test alpine.nonconvex_terms[[:(x[2]), :(x[2])]][:lifted_var_ref].args[2] == 6
        @test alpine.nonconvex_terms[[:(x[4]), :(x[4])]][:lifted_var_ref].args[2] == 5
        @test alpine.nonconvex_terms[[:(x[7]), :(x[7])]][:lifted_var_ref].args[2] == 8
    end

    @testset "Expression Parsing || Linear Lifting || general" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        )

        m = basic_linear_lift(solver = test_solver)

        alpine = _build(m) # Setup internal model

        lk = Vector{Any}(undef, 5)
        lk[1] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1), (-1.0, 2)])),
        )
        lk[2] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(3.0, 2), (-1.0, 3)])),
        )
        lk[3] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1), (1.0, 3)])),
        )
        lk[4] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1)])),
        )
        lk[5] = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 3.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2), (1.0, 1), (1.0, 3)])),
        )

        @test length(keys(alpine.linear_terms)) == 5
        yids = [8, 9, 12, 14, 16]
        for i in 1:length(keys(alpine.linear_terms))
            for j in keys(alpine.linear_terms)
                if alpine.linear_terms[j][:id] == i
                    @test j == lk[i]
                    @test j[:sign] == :+
                    @test j[:scalar] == lk[i][:scalar]
                    @test j[:coef_var] == lk[i][:coef_var]
                    @test alpine.linear_terms[j][:y_idx] == yids[i]
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

        @test alpine.nonconvex_terms[nlk1][:id] == 7
        @test alpine.nonconvex_terms[nlk2][:id] == 3
        @test alpine.nonconvex_terms[nlk3][:id] == 4
        @test alpine.nonconvex_terms[nlk4][:id] == 6
        @test alpine.nonconvex_terms[nlk5][:id] == 2
        @test alpine.nonconvex_terms[nlk6][:id] == 5
        @test alpine.nonconvex_terms[nlk7][:id] == 1
        @test alpine.nonconvex_terms[nlk8][:id] == 9
        @test alpine.nonconvex_terms[nlk9][:id] == 8

        @test alpine.nonconvex_terms[nlk1][:lifted_var_ref].args[2] == 13
        @test alpine.nonconvex_terms[nlk2][:lifted_var_ref].args[2] == 6
        @test alpine.nonconvex_terms[nlk3][:lifted_var_ref].args[2] == 7
        @test alpine.nonconvex_terms[nlk4][:lifted_var_ref].args[2] == 11
        @test alpine.nonconvex_terms[nlk5][:lifted_var_ref].args[2] == 5
        @test alpine.nonconvex_terms[nlk6][:lifted_var_ref].args[2] == 10
        @test alpine.nonconvex_terms[nlk7][:lifted_var_ref].args[2] == 4
        @test alpine.nonconvex_terms[nlk8][:lifted_var_ref].args[2] == 17
        @test alpine.nonconvex_terms[nlk9][:lifted_var_ref].args[2] == 15

        @test alpine.nonconvex_terms[nlk1][:nonlinear_type] == :MULTILINEAR
        @test alpine.nonconvex_terms[nlk2][:nonlinear_type] == :MONOMIAL
        @test alpine.nonconvex_terms[nlk3][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk4][:nonlinear_type] == :MONOMIAL
        @test alpine.nonconvex_terms[nlk5][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk6][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk7][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk8][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk9][:nonlinear_type] == :MONOMIAL

        @test alpine.nonconvex_terms[nlk1][:var_idxs] == [8, 9, 12]
        @test alpine.nonconvex_terms[nlk2][:var_idxs] == [2, 2]
        @test alpine.nonconvex_terms[nlk3][:var_idxs] == [2, 3]
        @test alpine.nonconvex_terms[nlk4][:var_idxs] == [8]
        @test alpine.nonconvex_terms[nlk5][:var_idxs] == [1, 3]
        @test alpine.nonconvex_terms[nlk6][:var_idxs] == [8, 9]
        @test alpine.nonconvex_terms[nlk7][:var_idxs] == [1, 2]
        @test alpine.nonconvex_terms[nlk8][:var_idxs] == [16, 15]
        @test alpine.nonconvex_terms[nlk9][:var_idxs] == [14]
    end

    @testset "Expression Parsing || Linear Lifting || brainpc3" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "disc_ratio" => 8,
            "log_level" => 100,
        )

        m = brainpc3(solver = test_solver)

        alpine = _build(m)

        @test alpine.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:y_idx] == 6913
        @test alpine.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:id] == 2
        @test alpine.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:lifted_constr_ref] ==
              :(x[6913] == (*)(x[6912]))
        @test alpine.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6912]), :(x[6912])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:y_idx] == 6970
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:id] == 25
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:lifted_constr_ref] ==
              :(x[6970] == x[6903] * x[6969])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6969])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:y_idx] == 6915
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:id] == 3
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:lifted_constr_ref] ==
              :(x[6915] == x[6903] * x[6914])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6914])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:y_idx] == 6920
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:id] == 5
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:lifted_constr_ref] ==
              :(x[6920] == x[6903] * x[6919])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6919])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:y_idx] == 6928
        @test alpine.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:id] == 8
        @test alpine.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:lifted_constr_ref] ==
              :(x[6928] == (*)(x[6927]))
        @test alpine.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6927]), :(x[6927])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:y_idx] == 6953
        @test alpine.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:id] == 18
        @test alpine.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:lifted_constr_ref] ==
              :(x[6953] == (*)(x[6952]))
        @test alpine.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6952]), :(x[6952])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:y_idx] == 6963
        @test alpine.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:id] == 22
        @test alpine.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:lifted_constr_ref] ==
              :(x[6963] == (*)(x[6962]))
        @test alpine.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6962]), :(x[6962])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:y_idx] == 6975
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:id] == 27
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:lifted_constr_ref] ==
              :(x[6975] == x[6903] * x[6974])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6974])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:y_idx] == 6935
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:id] == 11
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:lifted_constr_ref] ==
              :(x[6935] == x[6903] * x[6934])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6934])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:y_idx] == 6960
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:id] == 21
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:lifted_constr_ref] ==
              :(x[6960] == x[6903] * x[6959])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6959])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:y_idx] == 7015
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:id] == 43
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:lifted_constr_ref] ==
              :(x[7015] == x[6903] * x[7014])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7014])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:y_idx] == 6940
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:id] == 13
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:lifted_constr_ref] ==
              :(x[6940] == x[6903] * x[6939])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6939])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:y_idx] == 7018
        @test alpine.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:id] == 44
        @test alpine.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:lifted_constr_ref] ==
              :(x[7018] == (*)(x[7017]))
        @test alpine.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7017]), :(x[7017])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:y_idx] == 6930
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:id] == 9
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:lifted_constr_ref] ==
              :(x[6930] == x[6903] * x[6929])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6929])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:y_idx] == 7013
        @test alpine.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:id] == 42
        @test alpine.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:lifted_constr_ref] ==
              :(x[7013] == (*)(x[7012]))
        @test alpine.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7012]), :(x[7012])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:y_idx] == 6948
        @test alpine.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:id] == 16
        @test alpine.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:lifted_constr_ref] ==
              :(x[6948] == (*)(x[6947]))
        @test alpine.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6947]), :(x[6947])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:y_idx] == 7033
        @test alpine.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:id] == 50
        @test alpine.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:lifted_constr_ref] ==
              :(x[7033] == (*)(x[7032]))
        @test alpine.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7032]), :(x[7032])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:y_idx] == 6998
        @test alpine.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:id] == 36
        @test alpine.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:lifted_constr_ref] ==
              :(x[6998] == (*)(x[6997]))
        @test alpine.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6997]), :(x[6997])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:y_idx] == 7028
        @test alpine.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:id] == 48
        @test alpine.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:lifted_constr_ref] ==
              :(x[7028] == (*)(x[7027]))
        @test alpine.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7027]), :(x[7027])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:y_idx] == 7038
        @test alpine.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:id] == 52
        @test alpine.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:lifted_constr_ref] ==
              :(x[7038] == (*)(x[7037]))
        @test alpine.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7037]), :(x[7037])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:y_idx] == 7008
        @test alpine.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:id] == 40
        @test alpine.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:lifted_constr_ref] ==
              :(x[7008] == (*)(x[7007]))
        @test alpine.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7007]), :(x[7007])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:y_idx] == 6965
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:id] == 23
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:lifted_constr_ref] ==
              :(x[6965] == x[6903] * x[6964])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6964])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:y_idx] == 6925
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:id] == 7
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:lifted_constr_ref] ==
              :(x[6925] == x[6903] * x[6924])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6924])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:y_idx] == 6980
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:id] == 29
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:lifted_constr_ref] ==
              :(x[6980] == x[6903] * x[6979])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6979])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:y_idx] == 6995
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:id] == 35
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:lifted_constr_ref] ==
              :(x[6995] == x[6903] * x[6994])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6994])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:y_idx] == 6973
        @test alpine.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:id] == 26
        @test alpine.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:lifted_constr_ref] ==
              :(x[6973] == (*)(x[6972]))
        @test alpine.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6972]), :(x[6972])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:y_idx] == 7020
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:id] == 45
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:lifted_constr_ref] ==
              :(x[7020] == x[6903] * x[7019])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7019])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:y_idx] == 6993
        @test alpine.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:id] == 34
        @test alpine.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:lifted_constr_ref] ==
              :(x[6993] == (*)(x[6992]))
        @test alpine.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6992]), :(x[6992])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:y_idx] == 7003
        @test alpine.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:id] == 38
        @test alpine.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:lifted_constr_ref] ==
              :(x[7003] == (*)(x[7002]))
        @test alpine.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7002]), :(x[7002])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:y_idx] == 7023
        @test alpine.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:id] == 46
        @test alpine.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:lifted_constr_ref] ==
              :(x[7023] == (*)(x[7022]))
        @test alpine.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[7022]), :(x[7022])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:y_idx] == 6909
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:id] == 1
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:lifted_constr_ref] ==
              :(x[6909] == x[6903] * x[6908])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6908])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:y_idx] == 6950
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:id] == 17
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:lifted_constr_ref] ==
              :(x[6950] == x[6903] * x[6949])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6949])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:y_idx] == 7005
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:id] == 39
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:lifted_constr_ref] ==
              :(x[7005] == x[6903] * x[7004])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7004])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:y_idx] == 6918
        @test alpine.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:id] == 4
        @test alpine.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:lifted_constr_ref] ==
              :(x[6918] == (*)(x[6917]))
        @test alpine.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6917]), :(x[6917])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:y_idx] == 6938
        @test alpine.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:id] == 12
        @test alpine.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:lifted_constr_ref] ==
              :(x[6938] == (*)(x[6937]))
        @test alpine.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6937]), :(x[6937])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:y_idx] == 6933
        @test alpine.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:id] == 10
        @test alpine.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:lifted_constr_ref] ==
              :(x[6933] == (*)(x[6932]))
        @test alpine.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6932]), :(x[6932])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:y_idx] == 6990
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:id] == 33
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:lifted_constr_ref] ==
              :(x[6990] == x[6903] * x[6989])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6989])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:y_idx] == 6978
        @test alpine.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:id] == 28
        @test alpine.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:lifted_constr_ref] ==
              :(x[6978] == (*)(x[6977]))
        @test alpine.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6977]), :(x[6977])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:y_idx] == 6988
        @test alpine.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:id] == 32
        @test alpine.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:lifted_constr_ref] ==
              :(x[6988] == (*)(x[6987]))
        @test alpine.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6987]), :(x[6987])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:y_idx] == 7035
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:id] == 51
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:lifted_constr_ref] ==
              :(x[7035] == x[6903] * x[7034])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7034])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:y_idx] == 6968
        @test alpine.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:id] == 24
        @test alpine.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:lifted_constr_ref] ==
              :(x[6968] == (*)(x[6967]))
        @test alpine.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6967]), :(x[6967])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:y_idx] == 6958
        @test alpine.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:id] == 20
        @test alpine.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:lifted_constr_ref] ==
              :(x[6958] == (*)(x[6957]))
        @test alpine.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6957]), :(x[6957])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:y_idx] == 7010
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:id] == 41
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:lifted_constr_ref] ==
              :(x[7010] == x[6903] * x[7009])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7009])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:y_idx] == 6943
        @test alpine.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:id] == 14
        @test alpine.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:lifted_constr_ref] ==
              :(x[6943] == (*)(x[6942]))
        @test alpine.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6942]), :(x[6942])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:y_idx] == 6983
        @test alpine.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:id] == 30
        @test alpine.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:lifted_constr_ref] ==
              :(x[6983] == (*)(x[6982]))
        @test alpine.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6982]), :(x[6982])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:y_idx] == 7000
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:id] == 37
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:lifted_constr_ref] ==
              :(x[7000] == x[6903] * x[6999])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6999])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:y_idx] == 6945
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:id] == 15
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:lifted_constr_ref] ==
              :(x[6945] == x[6903] * x[6944])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6944])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:y_idx] == 7025
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:id] == 47
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:lifted_constr_ref] ==
              :(x[7025] == x[6903] * x[7024])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7024])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:y_idx] == 6955
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:id] == 19
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:lifted_constr_ref] ==
              :(x[6955] == x[6903] * x[6954])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6954])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:y_idx] == 6923
        @test alpine.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:id] == 6
        @test alpine.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:lifted_constr_ref] ==
              :(x[6923] == (*)(x[6922]))
        @test alpine.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:nonlinear_type] ==
              :MONOMIAL
        @test alpine.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6922]), :(x[6922])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:y_idx] == 6985
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:id] == 31
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:lifted_constr_ref] ==
              :(x[6985] == x[6903] * x[6984])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[6984])]][:constr_id] ==
              Set(Any[0])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:y_idx] == 7030
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:id] == 49
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:lifted_constr_ref] ==
              :(x[7030] == x[6903] * x[7029])
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:nonlinear_type] ==
              :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:y_type] == :Cont
        @test alpine.nonconvex_terms[Expr[:(x[6903]), :(x[7029])]][:constr_id] ==
              Set(Any[0])
        lk1 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6910), (-1.0, 6909)])),
        )
        @test alpine.linear_terms[lk1][:y_idx] == 6911
        @test alpine.linear_terms[lk1][:id] == 3
        @test alpine.linear_terms[lk1][:y_type] == :(Cont)
        lk2 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3467), (1.0, 16)])),
        )
        @test alpine.linear_terms[lk2][:y_idx] == 6908
        @test alpine.linear_terms[lk2][:id] == 1
        @test alpine.linear_terms[lk2][:y_type] == :(Cont)
        lk3 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(0.001548, 6903), (1.0, 7016), (1.0, 5702)]),
            ),
        )
        @test alpine.linear_terms[lk3][:y_idx] == 7017
        @test alpine.linear_terms[lk3][:id] == 67
        @test alpine.linear_terms[lk3][:y_type] == :(Cont)
        lk4 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 5162), (1.0, 1711)])),
        )
        @test alpine.linear_terms[lk4][:y_idx] == 7004
        @test alpine.linear_terms[lk4][:id] == 59
        @test alpine.linear_terms[lk4][:y_type] == :(Cont)
        lk5 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 7021), (0.001435, 6903), (1.0, 6002)]),
            ),
        )
        @test alpine.linear_terms[lk5][:y_idx] == 7022
        @test alpine.linear_terms[lk5][:id] == 70
        @test alpine.linear_terms[lk5][:y_type] == :(Cont)
        lk6 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6935), (1.0, 166)])),
        )
        @test alpine.linear_terms[lk6][:y_idx] == 6936
        @test alpine.linear_terms[lk6][:id] == 18
        @test alpine.linear_terms[lk6][:y_type] == :(Cont)
        lk7 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 4262), (0.00247, 6903), (1.0, 6981)]),
            ),
        )
        @test alpine.linear_terms[lk7][:y_idx] == 6982
        @test alpine.linear_terms[lk7][:id] == 46
        @test alpine.linear_terms[lk7][:y_type] == :(Cont)
        lk8 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7005), (1.0, 1711)])),
        )
        @test alpine.linear_terms[lk8][:y_idx] == 7006
        @test alpine.linear_terms[lk8][:id] == 60
        @test alpine.linear_terms[lk8][:y_type] == :(Cont)
        lk9 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1951), (1.0, 5402)])),
        )
        @test alpine.linear_terms[lk9][:y_idx] == 7009
        @test alpine.linear_terms[lk9][:id] == 62
        @test alpine.linear_terms[lk9][:y_type] == :(Cont)
        lk10 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 226), (-1.0, 6945)])),
        )
        @test alpine.linear_terms[lk10][:y_idx] == 6946
        @test alpine.linear_terms[lk10][:id] == 24
        @test alpine.linear_terms[lk10][:y_type] == :(Cont)
        lk11 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7030), (1.0, 3151)])),
        )
        @test alpine.linear_terms[lk11][:y_idx] == 7031
        @test alpine.linear_terms[lk11][:id] == 75
        @test alpine.linear_terms[lk11][:y_type] == :(Cont)
        lk12 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6991), (1.0, 4622), (0.002066, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk12][:y_idx] == 6992
        @test alpine.linear_terms[lk12][:id] == 52
        @test alpine.linear_terms[lk12][:y_type] == :(Cont)
        lk13 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 7026), (1.0, 6302), (0.001323, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk13][:y_idx] == 7027
        @test alpine.linear_terms[lk13][:id] == 73
        @test alpine.linear_terms[lk13][:y_type] == :(Cont)
        lk14 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 46), (1.0, 3497)])),
        )
        @test alpine.linear_terms[lk14][:y_idx] == 6914
        @test alpine.linear_terms[lk14][:id] == 5
        @test alpine.linear_terms[lk14][:y_type] == :(Cont)
        lk15 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6960), (1.0, 391)])),
        )
        @test alpine.linear_terms[lk15][:y_idx] == 6961
        @test alpine.linear_terms[lk15][:id] == 33
        @test alpine.linear_terms[lk15][:y_type] == :(Cont)
        lk16 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 631), (-1.0, 6975)])),
        )
        @test alpine.linear_terms[lk16][:y_idx] == 6976
        @test alpine.linear_terms[lk16][:id] == 42
        @test alpine.linear_terms[lk16][:y_type] == :(Cont)
        lk17 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1351), (-1.0, 6995)])),
        )
        @test alpine.linear_terms[lk17][:y_idx] == 6996
        @test alpine.linear_terms[lk17][:id] == 54
        @test alpine.linear_terms[lk17][:y_type] == :(Cont)
        lk18 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3557), (1.0, 106)])),
        )
        @test alpine.linear_terms[lk18][:y_idx] == 6924
        @test alpine.linear_terms[lk18][:id] == 11
        @test alpine.linear_terms[lk18][:y_type] == :(Cont)
        lk19 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 3647), (1.0, 6941), (0.004795, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk19][:y_idx] == 6942
        @test alpine.linear_terms[lk19][:id] == 22
        @test alpine.linear_terms[lk19][:y_type] == :(Cont)
        lk20 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3722), (1.0, 271)])),
        )
        @test alpine.linear_terms[lk20][:y_idx] == 6949
        @test alpine.linear_terms[lk20][:id] == 26
        @test alpine.linear_terms[lk20][:y_type] == :(Cont)
        lk21 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3647), (1.0, 196)])),
        )
        @test alpine.linear_terms[lk21][:y_idx] == 6939
        @test alpine.linear_terms[lk21][:id] == 20
        @test alpine.linear_terms[lk21][:y_type] == :(Cont)
        lk22 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 3467), (1.0, 6911), (6.34e-6, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk22][:y_idx] == 6912
        @test alpine.linear_terms[lk22][:id] == 4
        @test alpine.linear_terms[lk22][:y_type] == :(Cont)
        lk23 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1531), (1.0, 4982)])),
        )
        @test alpine.linear_terms[lk23][:y_idx] == 6999
        @test alpine.linear_terms[lk23][:id] == 56
        @test alpine.linear_terms[lk23][:y_type] == :(Cont)
        lk24 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 106), (-1.0, 6925)])),
        )
        @test alpine.linear_terms[lk24][:y_idx] == 6926
        @test alpine.linear_terms[lk24][:id] == 12
        @test alpine.linear_terms[lk24][:y_type] == :(Cont)
        lk25 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2251), (-1.0, 7015)])),
        )
        @test alpine.linear_terms[lk25][:y_idx] == 7016
        @test alpine.linear_terms[lk25][:id] == 66
        @test alpine.linear_terms[lk25][:y_type] == :(Cont)
        lk26 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6936), (0.004985, 6903), (1.0, 3617)]),
            ),
        )
        @test alpine.linear_terms[lk26][:y_idx] == 6937
        @test alpine.linear_terms[lk26][:id] == 19
        @test alpine.linear_terms[lk26][:y_type] == :(Cont)
        lk27 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6915), (1.0, 46)])),
        )
        @test alpine.linear_terms[lk27][:y_idx] == 6916
        @test alpine.linear_terms[lk27][:id] == 6
        @test alpine.linear_terms[lk27][:y_type] == :(Cont)
        lk28 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 3557), (0.005968, 6903), (1.0, 6926)]),
            ),
        )
        @test alpine.linear_terms[lk28][:y_idx] == 6927
        @test alpine.linear_terms[lk28][:id] == 13
        @test alpine.linear_terms[lk28][:y_type] == :(Cont)
        lk29 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 271), (-1.0, 6950)])),
        )
        @test alpine.linear_terms[lk29][:y_idx] == 6951
        @test alpine.linear_terms[lk29][:id] == 27
        @test alpine.linear_terms[lk29][:y_type] == :(Cont)
        lk30 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(0.002971, 6903), (1.0, 6971), (1.0, 3962)]),
            ),
        )
        @test alpine.linear_terms[lk30][:y_idx] == 6972
        @test alpine.linear_terms[lk30][:id] == 40
        @test alpine.linear_terms[lk30][:y_type] == :(Cont)
        lk31 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 511), (-1.0, 6970)])),
        )
        @test alpine.linear_terms[lk31][:y_idx] == 6971
        @test alpine.linear_terms[lk31][:id] == 39
        @test alpine.linear_terms[lk31][:y_type] == :(Cont)
        lk32 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 7031), (1.0, 6602), (0.00121, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk32][:y_idx] == 7032
        @test alpine.linear_terms[lk32][:id] == 76
        @test alpine.linear_terms[lk32][:y_type] == :(Cont)
        lk33 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(0.00166, 6903), (1.0, 7011), (1.0, 5402)]),
            ),
        )
        @test alpine.linear_terms[lk33][:y_idx] == 7012
        @test alpine.linear_terms[lk33][:id] == 64
        @test alpine.linear_terms[lk33][:y_type] == :(Cont)
        lk34 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, -0.003214),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 16)])),
        )
        @test alpine.linear_terms[lk34][:y_idx] == 6910
        @test alpine.linear_terms[lk34][:id] == 2
        @test alpine.linear_terms[lk34][:y_type] == :(Cont)
        lk35 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7035), (1.0, 3451)])),
        )
        @test alpine.linear_terms[lk35][:y_idx] == 7036
        @test alpine.linear_terms[lk35][:id] == 78
        @test alpine.linear_terms[lk35][:y_type] == :(Cont)
        lk36 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2551), (1.0, 6002)])),
        )
        @test alpine.linear_terms[lk36][:y_idx] == 7019
        @test alpine.linear_terms[lk36][:id] == 68
        @test alpine.linear_terms[lk36][:y_type] == :(Cont)
        lk37 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3842), (1.0, 391)])),
        )
        @test alpine.linear_terms[lk37][:y_idx] == 6959
        @test alpine.linear_terms[lk37][:id] == 32
        @test alpine.linear_terms[lk37][:y_type] == :(Cont)
        lk38 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6966), (0.00307, 6903), (1.0, 3902)]),
            ),
        )
        @test alpine.linear_terms[lk38][:y_idx] == 6967
        @test alpine.linear_terms[lk38][:id] == 37
        @test alpine.linear_terms[lk38][:y_type] == :(Cont)
        lk39 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 3722), (0.003829, 6903), (1.0, 6951)]),
            ),
        )
        @test alpine.linear_terms[lk39][:y_idx] == 6952
        @test alpine.linear_terms[lk39][:id] == 28
        @test alpine.linear_terms[lk39][:y_type] == :(Cont)
        lk40 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 451), (1.0, 3902)])),
        )
        @test alpine.linear_terms[lk40][:y_idx] == 6964
        @test alpine.linear_terms[lk40][:id] == 35
        @test alpine.linear_terms[lk40][:y_type] == :(Cont)
        lk41 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 7006), (1.0, 5162), (0.001713, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk41][:y_idx] == 7007
        @test alpine.linear_terms[lk41][:id] == 61
        @test alpine.linear_terms[lk41][:y_type] == :(Cont)
        lk42 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 811), (-1.0, 6980)])),
        )
        @test alpine.linear_terms[lk42][:y_idx] == 6981
        @test alpine.linear_terms[lk42][:id] == 45
        @test alpine.linear_terms[lk42][:y_type] == :(Cont)
        lk43 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 7001), (1.0, 4982), (0.00173, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk43][:y_idx] == 7002
        @test alpine.linear_terms[lk43][:id] == 58
        @test alpine.linear_terms[lk43][:y_type] == :(Cont)
        lk44 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6920), (1.0, 76)])),
        )
        @test alpine.linear_terms[lk44][:y_idx] == 6921
        @test alpine.linear_terms[lk44][:id] == 9
        @test alpine.linear_terms[lk44][:y_type] == :(Cont)
        lk45 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4622), (1.0, 1171)])),
        )
        @test alpine.linear_terms[lk45][:y_idx] == 6989
        @test alpine.linear_terms[lk45][:id] == 50
        @test alpine.linear_terms[lk45][:y_type] == :(Cont)
        lk46 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 991), (1.0, 4442)])),
        )
        @test alpine.linear_terms[lk46][:y_idx] == 6984
        @test alpine.linear_terms[lk46][:id] == 47
        @test alpine.linear_terms[lk46][:y_type] == :(Cont)
        lk47 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6930), (1.0, 136)])),
        )
        @test alpine.linear_terms[lk47][:y_idx] == 6931
        @test alpine.linear_terms[lk47][:id] == 15
        @test alpine.linear_terms[lk47][:y_type] == :(Cont)
        lk48 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 4802), (0.001891, 6903), (1.0, 6996)]),
            ),
        )
        @test alpine.linear_terms[lk48][:y_idx] == 6997
        @test alpine.linear_terms[lk48][:id] == 55
        @test alpine.linear_terms[lk48][:y_type] == :(Cont)
        lk49 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 511), (1.0, 3962)])),
        )
        @test alpine.linear_terms[lk49][:y_idx] == 6969
        @test alpine.linear_terms[lk49][:id] == 38
        @test alpine.linear_terms[lk49][:y_type] == :(Cont)
        lk50 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 811), (1.0, 4262)])),
        )
        @test alpine.linear_terms[lk50][:y_idx] == 6979
        @test alpine.linear_terms[lk50][:id] == 44
        @test alpine.linear_terms[lk50][:y_type] == :(Cont)
        lk51 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(0.01551, 6903), (1.0, 6916), (1.0, 3497)]),
            ),
        )
        @test alpine.linear_terms[lk51][:y_idx] == 6917
        @test alpine.linear_terms[lk51][:id] == 7
        @test alpine.linear_terms[lk51][:y_type] == :(Cont)
        lk52 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 451), (-1.0, 6965)])),
        )
        @test alpine.linear_terms[lk52][:y_idx] == 6966
        @test alpine.linear_terms[lk52][:id] == 36
        @test alpine.linear_terms[lk52][:y_type] == :(Cont)
        lk53 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(0.005773, 6903), (1.0, 3587), (1.0, 6931)]),
            ),
        )
        @test alpine.linear_terms[lk53][:y_idx] == 6932
        @test alpine.linear_terms[lk53][:id] == 16
        @test alpine.linear_terms[lk53][:y_type] == :(Cont)
        lk54 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 631), (1.0, 4082)])),
        )
        @test alpine.linear_terms[lk54][:y_idx] == 6974
        @test alpine.linear_terms[lk54][:id] == 41
        @test alpine.linear_terms[lk54][:y_type] == :(Cont)
        lk55 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 991), (-1.0, 6985)])),
        )
        @test alpine.linear_terms[lk55][:y_idx] == 6986
        @test alpine.linear_terms[lk55][:id] == 48
        @test alpine.linear_terms[lk55][:y_type] == :(Cont)
        lk56 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 4802), (1.0, 1351)])),
        )
        @test alpine.linear_terms[lk56][:y_idx] == 6994
        @test alpine.linear_terms[lk56][:id] == 53
        @test alpine.linear_terms[lk56][:y_type] == :(Cont)
        lk57 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6946), (0.004409, 6903), (1.0, 3677)]),
            ),
        )
        @test alpine.linear_terms[lk57][:y_idx] == 6947
        @test alpine.linear_terms[lk57][:id] == 25
        @test alpine.linear_terms[lk57][:y_type] == :(Cont)
        lk58 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 2251), (1.0, 5702)])),
        )
        @test alpine.linear_terms[lk58][:y_idx] == 7014
        @test alpine.linear_terms[lk58][:id] == 65
        @test alpine.linear_terms[lk58][:y_type] == :(Cont)
        lk59 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7020), (1.0, 2551)])),
        )
        @test alpine.linear_terms[lk59][:y_idx] == 7021
        @test alpine.linear_terms[lk59][:id] == 69
        @test alpine.linear_terms[lk59][:y_type] == :(Cont)
        lk60 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6940), (1.0, 196)])),
        )
        @test alpine.linear_terms[lk60][:y_idx] == 6941
        @test alpine.linear_terms[lk60][:id] == 21
        @test alpine.linear_terms[lk60][:y_type] == :(Cont)
        lk61 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 3527), (1.0, 6921), (0.008919, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk61][:y_idx] == 6922
        @test alpine.linear_terms[lk61][:id] == 10
        @test alpine.linear_terms[lk61][:y_type] == :(Cont)
        lk62 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7025), (1.0, 2851)])),
        )
        @test alpine.linear_terms[lk62][:y_idx] == 7026
        @test alpine.linear_terms[lk62][:id] == 72
        @test alpine.linear_terms[lk62][:y_type] == :(Cont)
        lk63 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6986), (1.0, 4442), (0.002252, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk63][:y_idx] == 6987
        @test alpine.linear_terms[lk63][:id] == 49
        @test alpine.linear_terms[lk63][:y_type] == :(Cont)
        lk64 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(0.002765, 6903), (1.0, 6976), (1.0, 4082)]),
            ),
        )
        @test alpine.linear_terms[lk64][:y_idx] == 6977
        @test alpine.linear_terms[lk64][:id] == 43
        @test alpine.linear_terms[lk64][:y_type] == :(Cont)
        lk65 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6956), (0.003269, 6903), (1.0, 3782)]),
            ),
        )
        @test alpine.linear_terms[lk65][:y_idx] == 6957
        @test alpine.linear_terms[lk65][:id] == 31
        @test alpine.linear_terms[lk65][:y_type] == :(Cont)
        lk66 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3587), (1.0, 136)])),
        )
        @test alpine.linear_terms[lk66][:y_idx] == 6929
        @test alpine.linear_terms[lk66][:id] == 14
        @test alpine.linear_terms[lk66][:y_type] == :(Cont)
        lk67 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 166), (1.0, 3617)])),
        )
        @test alpine.linear_terms[lk67][:y_idx] == 6934
        @test alpine.linear_terms[lk67][:id] == 17
        @test alpine.linear_terms[lk67][:y_type] == :(Cont)
        lk68 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 7010), (1.0, 1951)])),
        )
        @test alpine.linear_terms[lk68][:y_idx] == 7011
        @test alpine.linear_terms[lk68][:id] == 63
        @test alpine.linear_terms[lk68][:y_type] == :(Cont)
        lk69 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 7036), (0.001098, 6903), (1.0, 6902)]),
            ),
        )
        @test alpine.linear_terms[lk69][:y_idx] == 7037
        @test alpine.linear_terms[lk69][:id] == 79
        @test alpine.linear_terms[lk69][:y_type] == :(Cont)
        lk70 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3527), (1.0, 76)])),
        )
        @test alpine.linear_terms[lk70][:y_idx] == 6919
        @test alpine.linear_terms[lk70][:id] == 8
        @test alpine.linear_terms[lk70][:y_type] == :(Cont)
        lk71 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(
                :coef_var,
                Set(Any[(1.0, 6961), (1.0, 3842), (0.00317, 6903)]),
            ),
        )
        @test alpine.linear_terms[lk71][:y_idx] == 6962
        @test alpine.linear_terms[lk71][:id] == 34
        @test alpine.linear_terms[lk71][:y_type] == :(Cont)
        lk72 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6990), (1.0, 1171)])),
        )
        @test alpine.linear_terms[lk72][:y_idx] == 6991
        @test alpine.linear_terms[lk72][:id] == 51
        @test alpine.linear_terms[lk72][:y_type] == :(Cont)
        lk73 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3451), (1.0, 6902)])),
        )
        @test alpine.linear_terms[lk73][:y_idx] == 7034
        @test alpine.linear_terms[lk73][:id] == 77
        @test alpine.linear_terms[lk73][:y_type] == :(Cont)
        lk74 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6302), (1.0, 2851)])),
        )
        @test alpine.linear_terms[lk74][:y_idx] == 7024
        @test alpine.linear_terms[lk74][:id] == 71
        @test alpine.linear_terms[lk74][:y_type] == :(Cont)
        lk75 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1531), (-1.0, 7000)])),
        )
        @test alpine.linear_terms[lk75][:y_idx] == 7001
        @test alpine.linear_terms[lk75][:id] == 57
        @test alpine.linear_terms[lk75][:y_type] == :(Cont)
        lk76 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3151), (1.0, 6602)])),
        )
        @test alpine.linear_terms[lk76][:y_idx] == 7029
        @test alpine.linear_terms[lk76][:id] == 74
        @test alpine.linear_terms[lk76][:y_type] == :(Cont)
        lk77 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 3782), (1.0, 331)])),
        )
        @test alpine.linear_terms[lk77][:y_idx] == 6954
        @test alpine.linear_terms[lk77][:id] == 29
        @test alpine.linear_terms[lk77][:y_type] == :(Cont)
        lk78 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(-1.0, 6955), (1.0, 331)])),
        )
        @test alpine.linear_terms[lk78][:y_idx] == 6956
        @test alpine.linear_terms[lk78][:id] == 30
        @test alpine.linear_terms[lk78][:y_type] == :(Cont)
        lk79 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 226), (1.0, 3677)])),
        )
        @test alpine.linear_terms[lk79][:y_idx] == 6944
        @test alpine.linear_terms[lk79][:id] == 23
        @test alpine.linear_terms[lk79][:y_type] == :(Cont)
    end
end

@testset "Expression Parsing || Basic Multiplication Operators (Machine Generated for diffs)" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => JUNIPER,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )

    m = operator_basic(solver = test_solver)
    JuMP.set_optimize_hook(m, MOI.Utilities.attach_optimizer)
    JuMP.optimize!(m)
    JuMP.set_optimize_hook(m, nothing)
    alpine = JuMP.backend(m).optimizer.model
    Alpine.load!(alpine)

    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:y_idx] == 31
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:id] == 27
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:lifted_constr_ref] ==
          :(x[31] == x[3] * x[4])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:y_idx] == 83
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:id] == 79
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:lifted_constr_ref] ==
          :(x[83] == x[4] * x[5] * x[76])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[76])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:y_idx] == 9
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:id] == 5
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:lifted_constr_ref] ==
          :(x[9] == (*)(x[2]))
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[2])]][:nonlinear_type] == :MONOMIAL
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:y_idx] == 19
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:id] == 15
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:lifted_constr_ref] ==
          :(x[19] == x[5] * x[7])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:y_idx] == 32
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:id] == 28
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:lifted_constr_ref] ==
          :(x[32] == x[2] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:y_idx] == 49
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:id] == 45
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:lifted_constr_ref] ==
          :(x[49] == x[1] * x[2] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:y_idx] == 50
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:id] == 46
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:lifted_constr_ref] ==
          :(x[50] == x[3] * x[5] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5]), :(x[9])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:y_idx] == 69
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:id] == 65
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:lifted_constr_ref] ==
          :(x[69] == x[4] * x[47])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[47])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:y_idx] == 28
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:id] == 24
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:lifted_constr_ref] ==
          :(x[28] == (*)(x[4]))
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[4])]][:nonlinear_type] == :MONOMIAL
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:y_idx] == 37
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:id] == 33
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:lifted_constr_ref] ==
          :(x[37] == x[4] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:y_idx] == 5
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:id] == 1
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:lifted_constr_ref] ==
          :(x[5] == (*)(x[1]))
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[1])]][:nonlinear_type] == :MONOMIAL
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:y_idx] == 96
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:id] == 92
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:lifted_constr_ref] ==
          :(x[96] == x[5] * x[95])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[95])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:y_idx] == 41
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:id] == 37
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:lifted_constr_ref] ==
          :(x[41] == x[2] * x[40])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[40])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:y_idx] == 26
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:id] == 22
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:lifted_constr_ref] ==
          :(x[26] == x[6] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:y_idx] == 84
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:id] == 80
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:lifted_constr_ref] ==
          :(x[84] == x[1] * x[78] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[78]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:y_idx] == 58
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:id] == 54
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:lifted_constr_ref] ==
          :(x[58] == x[6] * x[43])
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[43])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:y_idx] == 106
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:id] == 102
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:lifted_constr_ref] ==
          :(x[106] == x[4] * x[6] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[6]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:y_idx] == 90
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:id] == 86
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:lifted_constr_ref] ==
          :(x[90] == x[2] * x[4] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:y_idx] == 42
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:id] == 38
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:lifted_constr_ref] ==
          :(x[42] == x[1] * x[41])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[41])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:y_idx] == 38
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:id] == 34
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:lifted_constr_ref] ==
          :(x[38] == x[2] * x[37])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[37])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:y_idx] == 87
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:id] == 83
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:lifted_constr_ref] ==
          :(x[87] == x[5] * x[85])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[85])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:y_idx] == 66
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:id] == 62
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:lifted_constr_ref] ==
          :(x[66] == x[2] * x[3] * x[5] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:y_idx] == 53
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:id] == 49
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:lifted_constr_ref] ==
          :(x[53] == x[5] * x[9] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:y_idx] == 107
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:id] == 103
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:lifted_constr_ref] ==
          :(x[107] == x[3] * x[6] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:y_idx] == 17
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:id] == 13
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:lifted_constr_ref] ==
          :(x[17] == x[2] * x[5])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:y_idx] == 57
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:id] == 53
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:lifted_constr_ref] ==
          :(x[57] == x[6] * x[37])
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[37])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:y_idx] == 81
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:id] == 77
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:lifted_constr_ref] ==
          :(x[81] == x[5] * x[7] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[7]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:y_idx] == 39
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:id] == 35
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:lifted_constr_ref] ==
          :(x[39] == x[1] * x[38])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[38])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:y_idx] == 48
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:id] == 44
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:lifted_constr_ref] ==
          :(x[48] == x[1] * x[3] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[9])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:y_idx] == 24
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:id] == 20
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:lifted_constr_ref] ==
          :(x[24] == x[4] * x[23])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[23])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:y_idx] == 23
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:id] == 19
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:lifted_constr_ref] ==
          :(x[23] == x[3] * x[17])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[17])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:y_idx] == 51
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:id] == 47
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:lifted_constr_ref] ==
          :(x[51] == x[1] * x[9] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:y_idx] == 75
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:id] == 71
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:lifted_constr_ref] ==
          :(x[75] == x[4] * x[5] * x[7])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[5]), :(x[7])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:y_idx] == 95
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:id] == 91
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] ==
          :(x[95] == x[9] * x[10] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:y_idx] == 11
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:id] == 7
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:lifted_constr_ref] ==
          :(x[11] == x[9] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:y_idx] == 100
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:id] == 96
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:lifted_constr_ref] ==
          :(x[100] == x[1] * x[2] * x[37])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[37])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:y_idx] == 62
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:id] == 58
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:lifted_constr_ref] ==
          :(x[62] == x[2] * x[3] * x[4] * x[5])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4]), :(x[5])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:y_idx] == 65
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:id] == 61
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:lifted_constr_ref] ==
          :(x[65] == x[1] * x[4] * x[9] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[9]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:y_idx] == 80
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:id] == 76
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:lifted_constr_ref] ==
          :(x[80] == x[1] * x[7] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:y_idx] == 46
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:id] == 42
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:lifted_constr_ref] ==
          :(x[46] == x[1] * x[2] * x[3])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:y_idx] == 59
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:id] == 55
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:lifted_constr_ref] ==
          :(x[59] == x[17] * x[40])
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[40])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:y_idx] == 86
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:id] == 82
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:lifted_constr_ref] ==
          :(x[86] == x[1] * x[85])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[85])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:y_idx] == 97
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:id] == 93
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:lifted_constr_ref] ==
          :(x[97] == x[1] * x[2] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[31])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:y_idx] == 29
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:id] == 25
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:lifted_constr_ref] ==
          :(x[29] == x[13] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:y_idx] == 92
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:id] == 88
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:lifted_constr_ref] ==
          :(x[92] == x[2] * x[3] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:y_idx] == 20
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:id] == 16
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:lifted_constr_ref] ==
          :(x[20] == x[1] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:y_idx] == 109
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:id] == 105
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:lifted_constr_ref] ==
          :(x[109] == x[14] * x[10] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:y_idx] == 93
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:id] == 89
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:lifted_constr_ref] ==
          :(x[93] == x[1] * x[92])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[92])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:y_idx] == 85
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:id] == 81
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] ==
          :(x[85] == x[2] * x[3] * x[4])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:y_idx] == 16
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:id] == 12
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:lifted_constr_ref] ==
          :(x[16] == x[1] * x[11])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:y_idx] == 25
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:id] == 21
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:lifted_constr_ref] ==
          :(x[25] == x[4] * x[21])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[21])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:y_idx] == 73
    @test alpine.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:id] == 69
    @test alpine.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:lifted_constr_ref] ==
          :(x[73] == x[53] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[53]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:y_idx] == 34
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:id] == 30
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:lifted_constr_ref] ==
          :(x[34] == x[5] * x[32])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[32])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:y_idx] == 40
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:id] == 36
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:lifted_constr_ref] ==
          :(x[40] == x[3] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:y_idx] == 67
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:id] == 63
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:lifted_constr_ref] ==
          :(x[67] == x[5] * x[9] * x[10] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[10]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:y_idx] == 102
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:id] == 98
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:lifted_constr_ref] ==
          :(x[102] == x[5] * x[9] * x[43])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9]), :(x[43])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:y_idx] == 27
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:id] == 23
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:lifted_constr_ref] ==
          :(x[27] == x[4] * x[26])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[26])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:y_idx] == 56
    @test alpine.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:id] == 52
    @test alpine.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:lifted_constr_ref] ==
          :(x[56] == x[20] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[20]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:y_idx] == 63
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:id] == 59
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:lifted_constr_ref] ==
          :(x[63] == x[1] * x[3] * x[4] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[3]), :(x[4]), :(x[9])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:y_idx] == 15
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:id] == 11
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:lifted_constr_ref] ==
          :(x[15] == x[14] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:y_idx] == 77
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:id] == 73
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:lifted_constr_ref] ==
          :(x[77] == x[1] * x[4] * x[76])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[76])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:y_idx] == 89
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:id] == 85
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:lifted_constr_ref] ==
          :(x[89] == x[1] * x[88])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[88])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:y_idx] == 94
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:id] == 90
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:lifted_constr_ref] ==
          :(x[94] == x[5] * x[90])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[90])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:y_idx] == 108
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:id] == 104
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:lifted_constr_ref] ==
          :(x[108] == x[6] * x[10] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[10]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:y_idx] == 70
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:id] == 66
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:lifted_constr_ref] ==
          :(x[70] == x[4] * x[48])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[48])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:y_idx] == 54
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:id] == 50
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:lifted_constr_ref] ==
          :(x[54] == x[6] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:y_idx] == 99
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:id] == 95
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:lifted_constr_ref] ==
          :(x[99] == x[1] * x[9] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[9]), :(x[31])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:y_idx] == 36
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:id] == 32
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:lifted_constr_ref] ==
          :(x[36] == x[1] * x[35])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[35])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:y_idx] == 14
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:id] == 10
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:lifted_constr_ref] ==
          :(x[14] == x[5] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:y_idx] == 22
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:id] == 18
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:lifted_constr_ref] ==
          :(x[22] == x[4] * x[13])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[13])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:y_idx] == 64
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:id] == 60
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:lifted_constr_ref] ==
          :(x[64] == x[1] * x[2] * x[10] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[10]), :(x[28])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:y_idx] == 72
    @test alpine.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:id] == 68
    @test alpine.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:lifted_constr_ref] ==
          :(x[72] == x[46] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[46]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:y_idx] == 55
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:id] == 51
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:lifted_constr_ref] ==
          :(x[55] == x[17] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:y_idx] == 98
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:id] == 94
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:lifted_constr_ref] ==
          :(x[98] == x[2] * x[5] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[31])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:y_idx] == 43
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:id] == 39
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:lifted_constr_ref] ==
          :(x[43] == x[10] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:y_idx] == 78
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:id] == 74
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:lifted_constr_ref] ==
          :(x[78] == x[2] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:y_idx] == 52
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:id] == 48
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:lifted_constr_ref] ==
          :(x[52] == x[2] * x[5] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[5]), :(x[10])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:y_idx] == 91
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:id] == 87
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:lifted_constr_ref] ==
          :(x[91] == x[1] * x[90])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[90])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:y_idx] == 104
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:id] == 100
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:lifted_constr_ref] ==
          :(x[104] == x[3] * x[4] * x[17])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[17])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:y_idx] == 76
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:id] == 72
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:lifted_constr_ref] ==
          :(x[76] == x[3] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:y_idx] == 79
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:id] == 75
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:lifted_constr_ref] ==
          :(x[79] == x[1] * x[4] * x[78])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[78])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:y_idx] == 10
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:id] == 6
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:lifted_constr_ref] ==
          :(x[10] == (*)(x[3]))
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[3])]][:nonlinear_type] == :MONOMIAL
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:y_idx] == 71
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:id] == 67
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:lifted_constr_ref] ==
          :(x[71] == x[4] * x[49])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[49])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:y_idx] == 88
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:id] == 84
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:lifted_constr_ref] ==
          :(x[88] == x[3] * x[4] * x[9])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[9])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:y_idx] == 30
    @test alpine.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:id] == 26
    @test alpine.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:lifted_constr_ref] ==
          :(x[30] == x[15] * x[28])
    @test alpine.nonconvex_terms[Expr[:(x[15]), :(x[28])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:y_idx] == 68
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:id] == 64
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:lifted_constr_ref] ==
          :(x[68] == x[4] * x[46])
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[46])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:y_idx] == 12
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:id] == 8
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:lifted_constr_ref] ==
          :(x[12] == x[5] * x[11])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[11])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:y_idx] == 74
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:id] == 70
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:lifted_constr_ref] ==
          :(x[74] == x[1] * x[4] * x[7])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[7])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:y_idx] == 105
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:id] == 101
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:lifted_constr_ref] ==
          :(x[105] == x[3] * x[4] * x[20])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[20])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:y_idx] == 82
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:id] == 78
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:lifted_constr_ref] ==
          :(x[82] == x[1] * x[4] * x[11])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[4]), :(x[11])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:y_idx] == 44
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:id] == 40
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:lifted_constr_ref] ==
          :(x[44] == x[9] * x[43])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[43])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:y_idx] == 61
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:id] == 57
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:lifted_constr_ref] ==
          :(x[61] == x[1] * x[2] * x[3] * x[4])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3]), :(x[4])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:y_idx] == 18
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:id] == 14
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:lifted_constr_ref] ==
          :(x[18] == x[17] * x[10])
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:y_idx] == 8
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:id] == 4
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:lifted_constr_ref] ==
          :(x[8] == x[1] * x[7])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 6
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:id] == 2
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:lifted_constr_ref] ==
          :(x[6] == x[1] * x[2])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:y_idx] == 101
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:id] == 97
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:lifted_constr_ref] ==
          :(x[101] == x[1] * x[2] * x[40])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[40])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:y_idx] == 60
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:id] == 56
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:lifted_constr_ref] ==
          :(x[60] == x[17] * x[37])
    @test alpine.nonconvex_terms[Expr[:(x[17]), :(x[37])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:y_idx] == 103
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:id] == 99
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:lifted_constr_ref] ==
          :(x[103] == x[3] * x[4] * x[6])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[4]), :(x[6])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:y_idx] == 7
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:id] == 3
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:lifted_constr_ref] ==
          :(x[7] == x[2] * x[3])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:y_idx] == 47
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:id] == 43
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:lifted_constr_ref] ==
          :(x[47] == x[2] * x[3] * x[5])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[3]), :(x[5])]][:nonlinear_type] ==
          :MULTILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:y_idx] == 13
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:id] == 9
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:lifted_constr_ref] ==
          :(x[13] == x[3] * x[6])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:y_idx] == 35
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:id] == 31
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:lifted_constr_ref] ==
          :(x[35] == x[9] * x[31])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[31])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:y_idx] == 33
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:id] == 29
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:lifted_constr_ref] ==
          :(x[33] == x[1] * x[32])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[32])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:y_idx] == 21
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:id] == 17
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:lifted_constr_ref] ==
          :(x[21] == x[3] * x[20])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[20])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:y_idx] == 45
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:id] == 41
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:lifted_constr_ref] ==
          :(x[45] == x[5] * x[44])
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[44])]][:nonlinear_type] == :BILINEAR

    @test alpine.bounding_constr_mip[1][:rhs] == 1.0
    @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[5])]
    @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[1][:sense] == :(>=)
    @test alpine.bounding_constr_mip[1][:cnt] == 1
    @test alpine.bounding_constr_mip[2][:rhs] == 1.0
    @test alpine.bounding_constr_mip[2][:vars] == Any[:(x[6])]
    @test alpine.bounding_constr_mip[2][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[2][:sense] == :(<=)
    @test alpine.bounding_constr_mip[2][:cnt] == 1
    @test alpine.bounding_constr_mip[3][:rhs] == 1.0
    @test alpine.bounding_constr_mip[3][:vars] == Any[:(x[5]), :(x[7])]
    @test alpine.bounding_constr_mip[3][:coefs] == Any[1.0, 1.0]
    @test alpine.bounding_constr_mip[3][:sense] == :(<=)
    @test alpine.bounding_constr_mip[3][:cnt] == 2
    @test alpine.bounding_constr_mip[4][:rhs] == 1.0
    @test alpine.bounding_constr_mip[4][:vars] == Any[:(x[8])]
    @test alpine.bounding_constr_mip[4][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[4][:sense] == :(>=)
    @test alpine.bounding_constr_mip[4][:cnt] == 1
    @test alpine.bounding_constr_mip[5][:rhs] == 1.0
    @test alpine.bounding_constr_mip[5][:vars] == Any[:(x[12])]
    @test alpine.bounding_constr_mip[5][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[5][:sense] == :(<=)
    @test alpine.bounding_constr_mip[5][:cnt] == 1
    @test alpine.bounding_constr_mip[6][:rhs] == 1.0
    @test alpine.bounding_constr_mip[6][:vars] == Any[:(x[13])]
    @test alpine.bounding_constr_mip[6][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[6][:sense] == :(>=)
    @test alpine.bounding_constr_mip[6][:cnt] == 1
    @test alpine.bounding_constr_mip[7][:rhs] == 1.0
    @test alpine.bounding_constr_mip[7][:vars] == Any[:(x[15])]
    @test alpine.bounding_constr_mip[7][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[7][:sense] == :(<=)
    @test alpine.bounding_constr_mip[7][:cnt] == 1
    @test alpine.bounding_constr_mip[8][:rhs] == 1.0
    @test alpine.bounding_constr_mip[8][:vars] == Any[:(x[16])]
    @test alpine.bounding_constr_mip[8][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[8][:sense] == :(>=)
    @test alpine.bounding_constr_mip[8][:cnt] == 1
    @test alpine.bounding_constr_mip[9][:rhs] == 1.0
    @test alpine.bounding_constr_mip[9][:vars] == Any[:(x[18])]
    @test alpine.bounding_constr_mip[9][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[9][:sense] == :(<=)
    @test alpine.bounding_constr_mip[9][:cnt] == 1
    @test alpine.bounding_constr_mip[10][:rhs] == 1.0
    @test alpine.bounding_constr_mip[10][:vars] == Any[:(x[19])]
    @test alpine.bounding_constr_mip[10][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[10][:sense] == :(>=)
    @test alpine.bounding_constr_mip[10][:cnt] == 1
    @test alpine.bounding_constr_mip[11][:rhs] == 1.0
    @test alpine.bounding_constr_mip[11][:vars] == Any[:(x[21])]
    @test alpine.bounding_constr_mip[11][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[11][:sense] == :(<=)
    @test alpine.bounding_constr_mip[11][:cnt] == 1
    @test alpine.bounding_constr_mip[12][:rhs] == 1.0
    @test alpine.bounding_constr_mip[12][:vars] == Any[:(x[22])]
    @test alpine.bounding_constr_mip[12][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[12][:sense] == :(>=)
    @test alpine.bounding_constr_mip[12][:cnt] == 1
    @test alpine.bounding_constr_mip[13][:rhs] == 1.0
    @test alpine.bounding_constr_mip[13][:vars] == Any[:(x[24])]
    @test alpine.bounding_constr_mip[13][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[13][:sense] == :(<=)
    @test alpine.bounding_constr_mip[13][:cnt] == 1
    @test alpine.bounding_constr_mip[14][:rhs] == 1.0
    @test alpine.bounding_constr_mip[14][:vars] == Any[:(x[25])]
    @test alpine.bounding_constr_mip[14][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[14][:sense] == :(>=)
    @test alpine.bounding_constr_mip[14][:cnt] == 1
    @test alpine.bounding_constr_mip[15][:rhs] == 1.0
    @test alpine.bounding_constr_mip[15][:vars] == Any[:(x[27])]
    @test alpine.bounding_constr_mip[15][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[15][:sense] == :(<=)
    @test alpine.bounding_constr_mip[15][:cnt] == 1
end

@testset "Expression Parsing || corner cases" begin
    @testset "Corner Cases - 1 : sign convertor special case" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "minlp_solver" => JUNIPER,
            "log_level" => 100,
        )

        m = Model(test_solver)
        @variable(m, 0 <= x[1:5] <= 1, Bin)
        @NLconstraint(m, x[1] + -x[2] >= 2)
        @objective(m, Min, x[1] + x[2])
        alpine = _build(m)

        @test alpine.bounding_constr_mip[1][:rhs] == 2.0
        @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0, -1.0]
        @test alpine.bounding_constr_mip[1][:sense] == :(>=)
        @test alpine.bounding_constr_mip[1][:cnt] == 2
    end

    @testset "Corner Cases - 2 : full sub-expression" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        )

        m = Model(test_solver)
        @variable(m, x[1:5] >= 0)
        @NLconstraint(m, (x[1] + 10 + 4 * x[2] + 200) * 5 >= 2 * 5)
        @constraint(m, x[1] + x[2] + 200 >= 2)
        @objective(m, Min, x[1] + x[2])
        alpine = _build(m)

        @test alpine.bounding_constr_mip[1][:rhs] == -198.0
        @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0]
        @test alpine.bounding_constr_mip[1][:sense] == :(>=)
        @test alpine.bounding_constr_mip[1][:cnt] == 2
        @test alpine.bounding_constr_mip[2][:rhs] == -1040.0
        @test alpine.bounding_constr_mip[2][:vars] == Any[:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[2][:coefs] == Any[5.0, 20.0]
        @test alpine.bounding_constr_mip[2][:sense] == :(>=)
        @test alpine.bounding_constr_mip[2][:cnt] == 2
    end

    @testset "Corner Cases - 2 : full sub-expression" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        )

        m = Model(test_solver)
        @variable(m, x[1:5] >= 0)
        @NLconstraint(m, (x[1] + 10 + 4 * x[2] + 200) * 5 >= 2 * 5)
        @constraint(m, x[1] + x[2] + 200 >= 2)
        @objective(m, Min, x[1] + x[2])
        alpine = _build(m)

        @test alpine.bounding_constr_mip[1][:rhs] == -198.0
        @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[1][:coefs] == Any[1.0, 1.0]
        @test alpine.bounding_constr_mip[1][:sense] == :(>=)
        @test alpine.bounding_constr_mip[1][:cnt] == 2
        @test alpine.bounding_constr_mip[2][:rhs] == -1040.0
        @test alpine.bounding_constr_mip[2][:vars] == Any[:(x[1]), :(x[2])]
        @test alpine.bounding_constr_mip[2][:coefs] == Any[5.0, 20.0]
        @test alpine.bounding_constr_mip[2][:sense] == :(>=)
        @test alpine.bounding_constr_mip[2][:cnt] == 2
    end
end

@testset "Expression Parsing || Discrete Multilinear" begin
    @testset "Expression Parsing || bmpl && binlin && binprod" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "minlp_solver" => JUNIPER,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        )

        m = bpml(solver = test_solver)

        alpine = _build(m) # Setup internal model

        @test length(keys(alpine.nonconvex_terms)) == 12

        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:y_idx] == 11
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:y_idx] == 12
        @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[6])]][:y_idx] == 13
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:y_idx] == 14
        @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[6])]][:y_idx] == 15
        @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:y_idx] == 16
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[16])]][:y_idx] == 17
        @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[8])]][:y_idx] == 18
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[18])]][:y_idx] == 19
        @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[18])]][:y_idx] == 20
        @test alpine.nonconvex_terms[Expr[
            :(x[1]),
            :(x[2]),
            :(x[3]),
            :(x[4]),
            :(x[5]),
        ]][:y_idx] == 21
        @test alpine.nonconvex_terms[Expr[
            :(x[6]),
            :(x[7]),
            :(x[9]),
            :(x[10]),
            :(x[6]),
        ]][:y_idx] == 22

        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2])]][:nonlinear_type] == :BINPROD
        @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[6])]][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[2]), :(x[3])]][:nonlinear_type] ==
              :BINPROD
        @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[6])]][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[16])]][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[7]), :(x[8])]][:nonlinear_type] ==
              :MULTILINEAR
        @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[18])]][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[Expr[:(x[14]), :(x[18])]][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[Expr[
            :(x[1]),
            :(x[2]),
            :(x[3]),
            :(x[4]),
            :(x[5]),
        ]][:nonlinear_type] == :BINPROD
        @test alpine.nonconvex_terms[Expr[
            :(x[6]),
            :(x[7]),
            :(x[9]),
            :(x[10]),
            :(x[6]),
        ]][:nonlinear_type] == :MULTILINEAR

        @test length(alpine.var_type) == 22
        @test alpine.var_type[11] == :Cont
        @test alpine.var_type[12] == :Bin
        @test alpine.var_type[13] == :Cont
        @test alpine.var_type[14] == :Bin
        @test alpine.var_type[15] == :Cont
        @test alpine.var_type[16] == :Cont
        @test alpine.var_type[17] == :Cont
        @test alpine.var_type[18] == :Cont
        @test alpine.var_type[19] == :Cont
        @test alpine.var_type[20] == :Cont
        @test alpine.var_type[21] == :Bin
        @test alpine.var_type[22] == :Cont
    end

    @testset "Expression Parsing || bmpl && binlin && binprod with linear lifting and coefficients" begin
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "minlp_solver" => JUNIPER,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "log_level" => 100,
        )

        m = bmpl_linearlifting(solver = test_solver)

        alpine = _build(m)

        lk1 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 1), (1.0, 17)])),
        )
        lk2 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 11), (1.0, 12)])),
        )
        lk3 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 7), (2.0, 6)])),
        )
        lk4 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6), (1.0, 19)])),
        )
        lk5 = Dict{Symbol,Any}(
            Pair{Symbol,Any}(:sign, :+),
            Pair{Symbol,Any}(:scalar, 0.0),
            Pair{Symbol,Any}(:coef_var, Set(Any[(1.0, 6), (1.0, 25)])),
        )

        @test haskey(alpine.linear_terms, lk1)
        @test haskey(alpine.linear_terms, lk2)
        @test haskey(alpine.linear_terms, lk3)
        @test haskey(alpine.linear_terms, lk4)
        @test haskey(alpine.linear_terms, lk5)

        @test alpine.linear_terms[lk1][:lifted_constr_ref] == :(x[20] == x[1] + x[17])
        @test alpine.linear_terms[lk2][:lifted_constr_ref] == :(x[13] == x[11] + x[12])
        @test alpine.linear_terms[lk3][:lifted_constr_ref] == :(x[15] == 2 * x[6] + x[7])
        @test alpine.linear_terms[lk4][:lifted_constr_ref] == :(x[21] == x[6] + x[19])
        @test alpine.linear_terms[lk5][:lifted_constr_ref] == :(x[26] == x[6] + x[25])

        @test alpine.linear_terms[lk1][:y_idx] == 20
        @test alpine.linear_terms[lk2][:y_idx] == 13
        @test alpine.linear_terms[lk3][:y_idx] == 15
        @test alpine.linear_terms[lk4][:y_idx] == 21
        @test alpine.linear_terms[lk5][:y_idx] == 26

        @test alpine.linear_terms[lk1][:y_type] == :Cont
        @test alpine.linear_terms[lk2][:y_type] == :Cont
        @test alpine.linear_terms[lk3][:y_type] == :Cont
        @test alpine.linear_terms[lk4][:y_type] == :Cont
        @test alpine.linear_terms[lk5][:y_type] == :Cont

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

        @test haskey(alpine.nonconvex_terms, nlk1)
        @test haskey(alpine.nonconvex_terms, nlk2)
        @test haskey(alpine.nonconvex_terms, nlk3)
        @test haskey(alpine.nonconvex_terms, nlk4)
        @test haskey(alpine.nonconvex_terms, nlk5)
        @test haskey(alpine.nonconvex_terms, nlk6)
        @test haskey(alpine.nonconvex_terms, nlk7)
        @test haskey(alpine.nonconvex_terms, nlk8)
        @test haskey(alpine.nonconvex_terms, nlk9)
        @test haskey(alpine.nonconvex_terms, nlk10)
        @test haskey(alpine.nonconvex_terms, nlk11)
        @test haskey(alpine.nonconvex_terms, nlk12)

        @test alpine.nonconvex_terms[nlk1][:id] == 7
        @test alpine.nonconvex_terms[nlk2][:id] == 1
        @test alpine.nonconvex_terms[nlk3][:id] == 9
        @test alpine.nonconvex_terms[nlk4][:id] == 12
        @test alpine.nonconvex_terms[nlk5][:id] == 3
        @test alpine.nonconvex_terms[nlk6][:id] == 4
        @test alpine.nonconvex_terms[nlk7][:id] == 5
        @test alpine.nonconvex_terms[nlk8][:id] == 11
        @test alpine.nonconvex_terms[nlk9][:id] == 6
        @test alpine.nonconvex_terms[nlk10][:id] == 2
        @test alpine.nonconvex_terms[nlk11][:id] == 8
        @test alpine.nonconvex_terms[nlk12][:id] == 10

        @test alpine.nonconvex_terms[nlk1][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk2][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk3][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk4][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk5][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk6][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk7][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk8][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk9][:y_type] == :Bin
        @test alpine.nonconvex_terms[nlk10][:y_type] == :Bin
        @test alpine.nonconvex_terms[nlk11][:y_type] == :Cont
        @test alpine.nonconvex_terms[nlk12][:y_type] == :Bin

        @test alpine.nonconvex_terms[nlk1][:y_idx] == 19
        @test alpine.nonconvex_terms[nlk2][:y_idx] == 11
        @test alpine.nonconvex_terms[nlk3][:y_idx] == 23
        @test alpine.nonconvex_terms[nlk4][:y_idx] == 27
        @test alpine.nonconvex_terms[nlk5][:y_idx] == 14
        @test alpine.nonconvex_terms[nlk6][:y_idx] == 16
        @test alpine.nonconvex_terms[nlk7][:y_idx] == 17
        @test alpine.nonconvex_terms[nlk8][:y_idx] == 25
        @test alpine.nonconvex_terms[nlk9][:y_idx] == 18
        @test alpine.nonconvex_terms[nlk10][:y_idx] == 12
        @test alpine.nonconvex_terms[nlk11][:y_idx] == 22
        @test alpine.nonconvex_terms[nlk12][:y_idx] == 24

        @test alpine.nonconvex_terms[nlk1][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[nlk2][:nonlinear_type] == :MULTILINEAR
        @test alpine.nonconvex_terms[nlk3][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk4][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[nlk5][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[nlk6][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[nlk7][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[nlk8][:nonlinear_type] == :BINLIN
        @test alpine.nonconvex_terms[nlk9][:nonlinear_type] == :BINPROD
        @test alpine.nonconvex_terms[nlk10][:nonlinear_type] == :BINPROD
        @test alpine.nonconvex_terms[nlk11][:nonlinear_type] == :BILINEAR
        @test alpine.nonconvex_terms[nlk12][:nonlinear_type] == :BINPROD
    end
end

@testset "Expr dereferencing for @NLexpression" begin
    # Taken from Juniper.jl
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => JUNIPER,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => true,
    )

    m = Model(test_solver)
    @variable(m, 0 <= x <= 1)
    @variable(m, 0 <= z <= 1, Int)
    b_expr = @NLexpression(m, z / 1)
    @NLconstraint(m, 1.1 >= b_expr)
    an_expr = @NLexpression(m, x / 1)
    @NLobjective(m, Max, z + an_expr)
    JuMP.optimize!(m)

    @test isapprox(JuMP.objective_value(m), 2.0; atol = 1E-6)
    @test isapprox(JuMP.value(z), 1; atol = 1e-6)
    @test isapprox(JuMP.value(x), 1; atol = 1e-6)
end
