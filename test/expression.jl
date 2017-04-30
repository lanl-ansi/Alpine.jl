@testset "Expression Parser Test" begin

    @testset " Expression Test || Bilinear || Simple || bi1.jl " begin

        bi1 = pod_example_bi1()
        bi1_d = JuMP.NLPEvaluator(bi1)
        MathProgBase.initialize(bi1_d, [:ExprGraph])
        bi1_pod = POD._pod_formulate_liftmodel(bi1_d)

        solve(bi1)
        solve(bi1_pod)

        @test isapprox(bi1.objVal, 222.0744; atol = 1e-3)
        @test isapprox(bi1_pod.objVal, 222.0744; atol = 1e-3)
        @test isa(bi1_pod.ext[:map], Dict)

        nlp1=pod_example_nlp1()
        nlp1_d = JuMP.NLPEvaluator(nlp1)
        MathProgBase.initialize(nlp1_d, [:ExprGraph])
        nlp1_pod = POD._pod_formulate_liftmodel(nlp1_d)

        solve(nlp1)
        solve(nlp1_pod)

        @test isapprox(nlp1.objVal, 58.38366169; atol = 1e-3)
        @test isapprox(nlp1_pod.objVal, 58.38366169; atol = 1e-3)
        @test isa(nlp1_pod.ext[:map], Dict)

        nlp3=pod_example_nlp3()
        nlp3_d = JuMP.NLPEvaluator(nlp3)
        MathProgBase.initialize(nlp3_d, [:ExprGraph])
        nlp3_pod = POD._pod_formulate_liftmodel(nlp3_d)

        nlp3_status = solve(nlp3)
        nlp3_pod_status = solve(nlp3_pod)

        @test nlp3_status == :Infeasible
        @test nlp3_pod_status == :Infeasible
        @test isa(nlp3_pod.ext[:map], Dict)
    end
end
