@testset "Expression Parser Test" begin

    @testset "Expression Test || bilinear || Affine || exprs.jl" begin
        m=example_exprs()
        solve(m)
        @test isapprox(m.objVal, 306.8500; atol = 1e-3)

        # -1.0 * x[1] <= 109.0
        ex = m.internalModel.lifted_constr_expr_mip[1]
        coef, var, rhs = POD.expr_to_affine(ex)
        @show coef, var, rhs
        @test coef == [-1.0]
        @test var == [:(x[1])]
        @test isapprox(rhs, 109.0; atol = 1e-3)

        # (-1.0 * x[12] + 20.0 * x[13]) - 222.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[2]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [3.0,3.0,3.0,3.0]
        @test var == [:(x[8]),:(x[9]),:(x[10]),:(x[11])]
        @test isapprox(rhs, 111.0; atol = 1e-3)

        # -1.0 * x[12] - 115.0 <= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[3]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [-1.0,20.0]
        @test var == [:(x[12]),:(x[13])]
        @test isapprox(rhs, 222.0; atol = 1e-3)

        # 1.0 * x[12] - 115.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[4]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [-1.0]
        @test var == [:(x[12])]     
        @test isapprox(rhs, 115.0; atol = 1e-3)

        # 1.0 * x[12] - 115.0 <= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[5]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [1.0]
        @test var == [:(x[12])]
        @test isapprox(rhs, 115.0; atol = 1e-3)

        # -1.0 * x[12] - 115.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[6]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [-1.0]
        @test var == [:(x[12])]
        @test isapprox(rhs, 115.0; atol = 1e-3)

        # (x[1] + 1.0 * x[14]) - 555.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[7]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [1.0, 1.0]
        @test var == [:(x[1]),:(x[14])]
        @test isapprox(rhs, 555.0; atol = 1e-3)

        # ((x[8] - 7.0 * x[9]) + x[10] + x[4]) - 6666.0 <= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[8]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [1.0,-7.0,1.0,1.0]
        @test var == [:(x[8]),:(x[9]),:(x[10]),:(x[4])]
        @test isapprox(rhs, 6666.0; atol = 1e-3)

        # ((13.0 * x[1] - x[2]) + 30.0 * x[3] + x[4]) - 77.0 >= 0.0
        ex = m.internalModel.lifted_constr_expr_mip[9]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [13.0,-1.0,30.0,1.0]
        @test var == [:(x[1]),:(x[2]),:(x[3]),:(x[4])]
        @test isapprox(rhs, 77.0; atol = 1e-3)
    
    end

    @testset "Expression Test || bilinear || Affine || nlp1.jl" begin

        m=example_nlp1()
        solve(m)

        ex = m.internalModel.lifted_constr_expr_mip[1]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [1.0]
        @test var == [:(x[5])]
        @test isapprox(rhs, 8.0; atol = 1e-3)
    
    end

    @testset "Expression Test || bilinear || Affine || nlp3.jl" begin

        m=example_nlp3()
        solve(m)

        ex = m.internalModel.lifted_constr_expr_mip[1]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [0.0025,0.0025]
        @test var == [:(x[4]),:(x[6])]
        @test isapprox(rhs, 1.0; atol = 1e-3)

        ex = m.internalModel.lifted_constr_expr_mip[2]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [-0.0025,0.0025,0.0025]
        @test var == [:(x[4]),:(x[5]),:(x[7])]
        @test isapprox(rhs, 1.0; atol = 1e-3)

        ex = m.internalModel.lifted_constr_expr_mip[3]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [-0.01, 0.01]
        @test var == [:(x[5]),:(x[8])]
        @test isapprox(rhs, 0.0; atol = 1e-3)

        ex = m.internalModel.lifted_constr_expr_mip[4]
        coef, var, rhs = POD.expr_to_affine(ex)
        @show coef
        @test (coef .== Any[100.0, -1.0, 833.33252]) == [true, true, true]
        @test var == [:(x[1]),:(x[9]),:(x[4])]
        @test isapprox(rhs, 83333.333; atol = 1e-3)

        ex = m.internalModel.lifted_constr_expr_mip[5]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [1.0,-1.0,-1250.0,1250.0]
        @test var == [:(x[10]),:(x[11]),:(x[4]),:(x[5])]
        @test isapprox(rhs, 0.0; atol = 1e-3)

        ex = m.internalModel.lifted_constr_expr_mip[6]
        coef, var, rhs = POD.expr_to_affine(ex)
        @test coef == [1.0,-1.0,-2500.0]
        @test var == [:(x[12]),:(x[13]),:(x[5])]
        @test isapprox(rhs, -1.25e6; atol = 1e-3)
        
    end

    @testset "Expression Test || bilinear || Simple || bi1.jl " begin

        bi1 = example_bi1()
        solve(bi1) # Setup internal model
        POD.populate_dict_nonlinear_info(bi1.internalModel)

        @test bi1.internalModel.dict_nonlinear_info[:xdim] == 6
        @test length(keys(bi1.internalModel.dict_nonlinear_info)) == 9
        @test haskey(bi1.internalModel.dict_nonlinear_info, [Expr(:ref, :x, 1), Expr(:ref, :x, 1)]) 
        @test haskey(bi1.internalModel.dict_nonlinear_info, [Expr(:ref, :x, 2), Expr(:ref, :x, 2)]) 
        @test haskey(bi1.internalModel.dict_nonlinear_info, [Expr(:ref, :x, 3), Expr(:ref, :x, 3)]) 
        @test haskey(bi1.internalModel.dict_nonlinear_info, [Expr(:ref, :x, 4), Expr(:ref, :x, 4)])
        @test haskey(bi1.internalModel.dict_nonlinear_info, [Expr(:ref, :x, 2), Expr(:ref, :x, 3)])
        @test haskey(bi1.internalModel.dict_nonlinear_info, [Expr(:ref, :x, 3), Expr(:ref, :x, 4)]) 

    end
end
