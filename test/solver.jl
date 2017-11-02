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
        test_solver = PODSolver(minlp_local_solver=BonminNLSolver(["bonmin.num_resolve_at_root=2";"bonmin.time_limit=60"]),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=0,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 13.359400664712709; atol=1e-3)
        @test isapprox(m.objBound, 13.977883333333327; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 26
        @test m.internalModel.discretization_var_pick_algo == 0

        # Minimum vertex cover
        test_solver = PODSolver(minlp_local_solver=BonminNLSolver(["bonmin.num_resolve_at_root=2";"bonmin.time_limit=60"]),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=1,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 13.359400664712709; atol=1e-3)
        @test isapprox(m.objBound, 14.0064; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 1

        # Adaptive Scheme vertex cover
        test_solver = PODSolver(minlp_local_solver=BonminNLSolver(["bonmin.num_resolve_at_root=2";"bonmin.time_limit=60"]),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=2,
                                discretization_uniform_rate=10,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 13.359400664712709; atol=1e-3)
        @test isapprox(m.objBound, 14.0064; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 2

        # Dynamic scheme
        test_solver = PODSolver(minlp_local_solver=BonminNLSolver(["bonmin.num_resolve_at_root=2";"bonmin.time_limit=60"]),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=3,
                                max_iter=1,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 13.359400664712709; atol=1e-3)
        @test isapprox(m.objBound, 14.0064; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.discretization_var_pick_algo == 3

        # Dynamic Scheme step 2
        test_solver = PODSolver(minlp_local_solver=BonminNLSolver(["bonmin.num_resolve_at_root=2";"bonmin.time_limit=60"]),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=3,
                                max_iter=2,
                                log_level=1)

        m = blend029_gl(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(m.objVal, 13.359400664712709; atol=1e-3)
        @test isapprox(m.objBound, 14.0064; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.var_discretization_mip == Any[47, 40, 46, 43, 60, 44, 37, 38, 57, 41]
        @test m.internalModel.discretization_var_pick_algo == 3

        # Dynamic Scheme step 2
        test_solver = PODSolver(minlp_local_solver=BonminNLSolver(["bonmin.num_resolve_at_root=2";"bonmin.time_limit=60"]),
                                nlp_local_solver=IpoptSolver(print_level=0),
                                mip_solver=CbcSolver(logLevel=0),
                                discretization_var_pick_algo=3,
                                log_level=100)

        m = blend029_gl(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 13.359400664712709; atol=1e-3)
        @test isapprox(m.objBound, 13.35940; atol=1e-3)

        @test length(m.internalModel.all_nonlinear_vars) == 26
        @test m.internalModel.all_nonlinear_vars == [26, 27, 29, 30, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60]
        @test length(m.internalModel.var_discretization_mip) == 10
        @test m.internalModel.var_discretization_mip == Any[47, 40, 46, 43, 60, 44, 37, 38, 57, 41]
        @test m.internalModel.discretization_var_pick_algo == 3
    end


end
