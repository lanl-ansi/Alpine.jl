@testset "Expression Parser Test" begin
    @testset " Expression Test || Bilinear || Simple || bi1.jl " begin

        include("../examples/bi1.jl")
        include("../src/nlexpr.jl")

        m = pod_example_bi1()
        d = JuMP.NLPEvaluator(m)
        MathProgBase.initialize(d, [:ExprGraph])

        objexpr = mpb.obj_expr(d)
        consexprs = []

        for i in 1:length(d.constraints)
            push!(consexprs, mpb.constr_expr(d,i).args[2])
        end

        A, B, C = pod_expr_rebuild(objexpr)
        for i in 1:length(d.constraints)
            A, B, C = pod_expr_rebuild(consexprs[i], A, B)
        end

        @test length(A) == 12
        @test length(B) == 8
        @test ([Expr(:ref, :x, 1)] in A)
        @test ([Expr(:ref, :x, 2)] in A)
        @test ([Expr(:ref, :x, 3)] in A)
        @test ([Expr(:ref, :x, 4)] in A)
        @test ([Expr(:ref, :x, 1), Expr(:ref, :x, 1)] in A)
        @test ([Expr(:ref, :x, 2), Expr(:ref, :x, 2)] in A)
        @test ([Expr(:ref, :x, 3), Expr(:ref, :x, 3)] in A)
        @test ([Expr(:ref, :x, 4), Expr(:ref, :x, 4)] in A)
        @test ([Expr(:ref, :x, 1), Expr(:ref, :x, 4)] in A) || ([Expr(:ref, :x, 4), Expr(:ref, :x, 1)] in A)
        @test ([Expr(:ref, :x, 1), Expr(:ref, :x, 2)] in A) || ([Expr(:ref, :x, 2), Expr(:ref, :x, 1)] in A)
        @test ([Expr(:ref, :x, 3), Expr(:ref, :x, 4)] in A) || ([Expr(:ref, :x, 4), Expr(:ref, :x, 3)] in A)
        @test ([Expr(:ref, :x, 2), Expr(:ref, :x, 3)] in A) || ([Expr(:ref, :x, 3), Expr(:ref, :x, 2)] in A)

        @test ([Expr(:ref, :x, 1), Expr(:ref, :x, 1)] in keys(B))
        @test ([Expr(:ref, :x, 2), Expr(:ref, :x, 2)] in keys(B))
        @test ([Expr(:ref, :x, 3), Expr(:ref, :x, 3)] in keys(B))
        @test ([Expr(:ref, :x, 4), Expr(:ref, :x, 4)] in keys(B))
        @test ([Expr(:ref, :x, 1), Expr(:ref, :x, 4)] in keys(B)) || ([Expr(:ref, :x, 4), Expr(:ref, :x, 1)] in keys(B))
        @test ([Expr(:ref, :x, 1), Expr(:ref, :x, 2)] in keys(B)) || ([Expr(:ref, :x, 2), Expr(:ref, :x, 1)] in keys(B))
        @test ([Expr(:ref, :x, 3), Expr(:ref, :x, 4)] in keys(B)) || ([Expr(:ref, :x, 4), Expr(:ref, :x, 3)] in keys(B))
        @test ([Expr(:ref, :x, 2), Expr(:ref, :x, 3)] in keys(B)) || ([Expr(:ref, :x, 3), Expr(:ref, :x, 2)] in keys(B))

        # More tests to be added
        # Perform a throughal test below

        pod_ys = []
        for i in keys(B)
            push!(pod_ys, B[i])
        end
        @test length(unique(pod_ys)) == 8

    end
end
