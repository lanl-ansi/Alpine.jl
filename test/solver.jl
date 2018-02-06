@testset "Solver Function Tests" begin

    @testset "PODNonlinearModel loading tests" begin
        # Random Model 1
        test_solver = PODSolver(nlp_solver=IpoptSolver(),
    						   mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = operator_c(solver=test_solver)

        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)

        # Expression Model 1
        test_solver = PODSolver(nlp_solver=IpoptSolver(),
    						   mip_solver=CbcSolver(logLevel=0),loglevel=100)
        m = exprstest(solver=test_solver)
        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)
    end

    @testset "Partitioning variable selection tests :: nlp3" begin

        # Select all NL variable
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=0,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound, 3004.2470074351413;atol=1e-3)
        @test length(m.internalModel.candidate_disc_vars) == 8
        @test length(m.internalModel.disc_vars) == 8
        @test m.internalModel.disc_var_pick == 0

        # Select all NL variable
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=2,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound, 3004.2470074351413;atol=1e-3)
        @test length(m.internalModel.candidate_disc_vars) == 8
        @test length(m.internalModel.disc_vars) == 8
        @test m.internalModel.disc_var_pick == 2

        # Minimum vertex cover algorithm
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=1,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound,2606.2285443624664;atol=1e-3)
        @test length(m.internalModel.candidate_disc_vars) == 8
        @test length(m.internalModel.disc_vars) == 3
        @test m.internalModel.disc_var_pick == 1

        # Adaptive variable selection scheme :: disc_var_pick = 3
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=3,
                                presolve_bp = false,
                                maxiter=2,
                                loglevel=100)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound, 4896.6075;atol=1e-3)
        @test length(m.internalModel.candidate_disc_vars) == 8
        @test length(m.internalModel.disc_vars) == 8
        @test m.internalModel.disc_var_pick == 3
    end

    @testset "Partitioning variable selection tests :: castro2m2" begin

        # Select all NL variable
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=0,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)

        m = castro2m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test m.objVal <= 470.3176
        @test isapprox(m.objBound, 77.9999999999999; atol=1e-3)

        @test length(m.internalModel.candidate_disc_vars) == 10
        @test length(m.internalModel.disc_vars) == 10
        @test m.internalModel.disc_var_pick == 0

        # Select minimum vertex cover
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=1,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)

        m = castro2m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test m.objVal <= 470.3176
        @test isapprox(m.objBound, 250055.0761; atol=1e-3)
        @test length(m.internalModel.candidate_disc_vars) == 10
        @test length(m.internalModel.disc_vars) == 4
        @test m.internalModel.disc_var_pick == 1

        # Criteria 15 static selection
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=2,
                                disc_uniform_rate=15,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)

        m = castro2m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test m.objVal <= 470.3176
        @test isapprox(m.objBound, 77.9999999999999; atol=1e-3)

        @test length(m.internalModel.candidate_disc_vars) == 10
        @test length(m.internalModel.disc_vars) == 10
        @test m.internalModel.disc_var_pick == 2
    end

    @testset "Partitioning variable selection tests :: blend029" begin

        # Select all NL variable
        test_solver = PODSolver(minlp_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                                nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=0,
                                disc_uniform_rate=10,
                                maxiter=1,
                                loglevel=100)

        m = blend029_gl(solver=test_solver)
        JuMP.build(m)

        @test length(m.internalModel.candidate_disc_vars) == 26
        @test Set(m.internalModel.candidate_disc_vars) == Set([26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60])
        @test length(m.internalModel.disc_vars) == 26
        @test m.internalModel.disc_var_pick == 0

        # Minimum vertex cover
        test_solver = PODSolver(minlp_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                                nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=1,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)

        m = blend029_gl(solver=test_solver)
        JuMP.build(m)

        @test length(m.internalModel.candidate_disc_vars) == 26
        @test Set(m.internalModel.candidate_disc_vars) == Set([26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60])
        @test length(m.internalModel.disc_vars) == 10
        @test m.internalModel.disc_var_pick == 1

        # Adaptive Scheme vertex cover
        test_solver = PODSolver(minlp_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                                nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=2,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)

        m = blend029_gl(solver=test_solver)
        JuMP.build(m)

        @test length(m.internalModel.candidate_disc_vars) == 26
        @test Set(m.internalModel.candidate_disc_vars) == Set([26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60])
        @test length(m.internalModel.disc_vars) == 10
        @test m.internalModel.disc_var_pick == 2
    end

    @testset "Partitioning variable selection tests :: castro6m2" begin

        # Dynamic Scheme step 2
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=3,
                                presolve_bp=true,
                                maxiter=1,
                                loglevel=100)

        m = castro6m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test m.objVal <= 228.7810
        @test isapprox(m.objBound, 106.05582679267336; atol=1e-3)

        @test length(m.internalModel.candidate_disc_vars) == 24
        @test Set(m.internalModel.candidate_disc_vars) == Set([26, 27, 28, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 122, 123, 124])
        @test length(m.internalModel.disc_vars) == 12
        @test Set(m.internalModel.disc_vars) == Set([122, 114, 109, 107, 123, 28, 110, 111, 112, 113, 108, 115])
        @test m.internalModel.disc_var_pick == 3

        # Dynamic Scheme step 2
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=3,
                                presolve_bp=true,
                                maxiter=2,
                                loglevel=100)

        m = castro6m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test m.objVal <= 228.7810
        @test m.objBound >= 127.2135

        @test length(m.internalModel.candidate_disc_vars) == 24
        @test Set(m.internalModel.candidate_disc_vars) == Set([26, 27, 28, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 122, 123, 124])
        @test length(m.internalModel.disc_vars) == 12
        @test Set(m.internalModel.disc_vars) == Set([124, 114, 101, 26, 102, 103, 123, 116, 118, 117, 113, 115])
        @test m.internalModel.disc_var_pick == 3
    end

    @testset "Test getsolvetime for time trackikng" begin
        test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                disc_var_pick=0,
                                disc_uniform_rate=10,
                                presolve_bp = false,
                                maxiter=1,
                                loglevel=100)

        m = castro2m2(solver=test_solver)
        status = solve(m)
        @test getsolvetime(m) > 0.
    end
end
