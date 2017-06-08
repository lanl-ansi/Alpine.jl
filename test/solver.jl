@testset "Solver Function Tests" begin

    using POD

    @testset "PODNonlinearModel loading tests" begin

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
