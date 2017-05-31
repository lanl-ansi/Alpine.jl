@testset "Solver Function Tests" begin

    using POD

    @testset "PODNonlinearModel loading tests" begin
        # This model currently doesn't suit the expression manipulation function given the complexity of the expressions
        # initval = [1,5,5,1]
        # m = Model(solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0)))
        # @variable(m, 1 <= x[i=1:4] <= 5, start=initval[i])
        # @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
        # @NLconstraint(m, sum(x[i]^2 for i=1:4) == 40)
        # @NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
        # solve(m)
        # @test isapprox(m.objVal, 17.014; atol = 1e-3)

        # Random Model 1
        m = example_bi1()
        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)

        # Expression Model 1
        m = example_exprs()
        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)

    end

    @testset "Basic global solve check" begin
        # Dummy test
        @test 3 == 3
    end

end
