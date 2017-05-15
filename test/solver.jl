@testset "Solver Function Tests" begin

    using POD

    @testset " Local Solve Tests => Function solve() " begin
        initval = [1,5,5,1]
        m = Model(solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0)))
        @variable(m, 1 <= x[i=1:4] <= 5, start=initval[i])
        @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
        @NLconstraint(m, sum(x[i]^2 for i=1:4) == 40)
        @NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
        solve(m)
        @test isapprox(m.objVal, 17.014; atol = 1e-3)

        # Random Model 1
        m = example_bi1()
        status = solve(m)
        @test status == :Optimal
        @test isapprox(m.objVal, 222.0744; atol = 1e-3)

        # NLP 1
        m = example_nlp1()
        status = solve(m)
        @test status == :Optimal
        @test isapprox(m.objVal, 58.38366169; atol = 1e-3)

        # NLP 2
        m = example_nlp2()
        status = solve(m)
        @test status == :Optimal
        @test isapprox(m.objVal, 5.00; atol = 1e-3)

        # NLP3
        m = example_nlp3()
        status = solve(m)
        @test status == :Infeasible
    end

end
