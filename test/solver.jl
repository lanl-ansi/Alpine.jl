@testset "Solver Function Tests" begin

    @testset "Dummy test" begin
        @test 3 == 3
    end

    @testset "PODNonlinearModel loading tests" begin
        # Random Model 1
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(),
    						   mip_solver=CbcSolver(logLevel=0),log_level=0)
        m = operator_c(solver=test_solver)

        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)

        # Expression Model 1
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(),
    						   mip_solver=CbcSolver(logLevel=0),log_level=0)
        m = exprstest(solver=test_solver)
        status = JuMP.build(m)
        @test isa(m.internalModel, POD.PODNonlinearModel)
    end

    @testset "Partitioning variable selection tests :: nlp3" begin

        # Select all NL variable
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=0,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound, 3004.2470074351413;atol=1e-3)
        @test length(m.internalModel.all_nonlinear_vars) == 8
        @test length(m.internalModel.var_discretization_mip) == 8
        @test m.internalModel.discretization_var_pick_algo == 0

        # Select all NL variable
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=2,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound, 3004.2470074351413;atol=1e-3)
        @test length(m.internalModel.all_nonlinear_vars) == 8
        @test length(m.internalModel.var_discretization_mip) == 8
        @test m.internalModel.discretization_var_pick_algo == 2

        # Minimum vertex cover algorithm
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=1,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound,2606.2285443624664;atol=1e-3)
        @test length(m.internalModel.all_nonlinear_vars) == 8
        @test length(m.internalModel.var_discretization_mip) == 3
        @test m.internalModel.discretization_var_pick_algo == 1

        # Adaptive variable selection scheme :: discretization_var_pick_algo = 3
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=3,
                                max_iter=2,
                                log_level=1)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 7049.2478976; atol=1e-3)
        @test isapprox(m.objBound, 4896.6075;atol=1e-3)
        @test length(m.internalModel.all_nonlinear_vars) == 8
        @test length(m.internalModel.var_discretization_mip) == 8
        @test m.internalModel.discretization_var_pick_algo == 3
    end

    @testset "Partitioning variable selection tests :: castro2m2" begin

        # Select all NL variable
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=0,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = castro2m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 130.70555147302088; atol=1e-3)
        @test isapprox(m.objBound, 77.9999999999999; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 10
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 0

        # Select minimum vertex cover
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=1,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = castro2m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 130.70555147302088; atol=1e-3)
        @test isapprox(m.objBound, 250055.0761; atol=1e-3)
        @test length(m.internalModel.all_nonlinear_vars) == 10
        @test length(m.internalModel.var_discretization_mip) == 4
        @test m.internalModel.discretization_var_pick_algo == 1

        # Criteria 15 static selection
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=2,
                                discretization_uniform_rate=15,
                                max_iter=1,
                                log_level=1)

        m = castro2m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 130.70555147302088; atol=1e-3)
        @test isapprox(m.objBound, 77.9999999999999; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 10
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 2
    end

    @testset "Partitioning variable selection tests :: blend029" begin

        # Select all NL variable
        test_solver = PODSolver(minlp_local_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=1),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=0,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        JuMP.build(m)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 26
        @test m.internalModel.discretization_var_pick_algo == 0

        # Minimum vertex cover
        test_solver = PODSolver(minlp_local_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=1),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=1,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        JuMP.build(m)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 1

        # Adaptive Scheme vertex cover
        test_solver = PODSolver(minlp_local_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=1),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=2,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        JuMP.build(m)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 2
    end

    @testset "Partitioning variable selection tests :: castro5m2" begin

        # Dynamic Scheme step 2
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(),
                                discretization_var_pick_algo=3,
                                bound_basic_propagation=true,
                                max_iter=1,
                                log_level=1)

        m = castro6m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 228.78085060535983; atol=1e-3)
        @test isapprox(m.objBound, 106.05582679267336; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 24
        @test m.internalModel.all_nonlinear_vars == [26, 27, 28, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 122, 123, 124]
        @test length(m.internalModel.var_discretization_mip) == 12
        @test m.internalModel.var_discretization_mip == Any[122, 114, 109, 107, 123, 28, 110, 111, 112, 113, 108, 115]
        @test m.internalModel.discretization_var_pick_algo == 3

        # Dynamic Scheme step 2
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=3,
                                bound_basic_propagation=true,
                                max_iter=2,
                                log_level=100)

        m = castro6m2(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 193.7271905260245; atol=1e-3)
        @test isapprox(m.objBound, 127.38069214825349; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 24
        @test m.internalModel.all_nonlinear_vars == [26, 27, 28, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 122, 123, 124]
        @test length(m.internalModel.var_discretization_mip) == 12
        @test m.internalModel.var_discretization_mip == Any[124, 114, 101, 26, 102, 103, 123, 116, 118, 117, 113, 115]
        @test m.internalModel.discretization_var_pick_algo == 3
    end

end
