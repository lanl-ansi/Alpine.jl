@testset "Expression Parser Test" begin

    @testset " Expression Test || Bilinear || Simple || bi1.jl " begin

        bi1 = pod_example_bi1()
        solve(bi1) # Setup internal model
        POD.populate_map(bi1.internalModel)

        @test bi1.internalModel.map_nonlinear_terms[:xdim] == 6
        @test length(keys(bi1.internalModel.map_nonlinear_terms)) == 8
        @test haskey(bi1.internalModel.map_nonlinear_terms, [Expr(:ref, :x, 1), Expr(:ref, :x, 1)]) 
        @test haskey(bi1.internalModel.map_nonlinear_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 2)]) 
        @test haskey(bi1.internalModel.map_nonlinear_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 3)]) 
        @test haskey(bi1.internalModel.map_nonlinear_terms, [Expr(:ref, :x, 4), Expr(:ref, :x, 4)])
        @test haskey(bi1.internalModel.map_nonlinear_terms, [Expr(:ref, :x, 2), Expr(:ref, :x, 3)])
        @test haskey(bi1.internalModel.map_nonlinear_terms, [Expr(:ref, :x, 3), Expr(:ref, :x, 4)]) 

    end
end
