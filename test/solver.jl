@testset "Solver Function Tests" begin

    using POD

    @testset "PODNonlinearModel loading tests" begin
        # Random Model 1
        m = operator_c()
        setsolver(m, PODSolver(nlp_local_solver=IpoptSolver(),
    						   mip_solver=CbcSolver(OutputFlag=0),log_level=0))
        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)

        # Expression Model 1
        m = exprstest()
        setsolver(m, PODSolver(nlp_local_solver=IpoptSolver(),
    						   mip_solver=CbcSolver(OutputFlag=0),log_level=0))
        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)

    end

    @testset "Basic global solve check" begin
        # Dummy test
        @test 3 == 3
    end

end
