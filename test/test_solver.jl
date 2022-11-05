@testset "Optimizer loading tests" begin
    # Random Model 1
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )
    m = operator_c(solver = test_solver)

    alpine = _build(m)
    @test isa(alpine, Alpine.Optimizer)

    # Expression Model 1
    m = exprstest(solver = test_solver)
    alpine = _build(m)
    @test isa(alpine, Alpine.Optimizer)
end

@testset "Partitioning variable selection tests :: nlp3" begin
    # Select all NL variable
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 0,
        "disc_uniform_rate" => 10,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "max_iter" => 1,
        "log_level" => 100,
    )

    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(JuMP.objective_value(m), 7049.2478976; atol = 1e-3)
    @test length(alpine.candidate_disc_vars) == 8
    @test length(alpine.disc_vars) == 8
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 0

    # Select all NL variable
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 2,
        "disc_uniform_rate" => 10,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "max_iter" => 1,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(JuMP.objective_value(m), 7049.2478976; atol = 1e-3)
    @test length(alpine.candidate_disc_vars) == 8
    @test length(alpine.disc_vars) == 8
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 2

    # Minimum vertex cover algorithm
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 1,
        "disc_uniform_rate" => 10,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "max_iter" => 1,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(JuMP.objective_value(m), 7049.2478976; atol = 1e-3)
    @test length(alpine.candidate_disc_vars) == 8
    @test length(alpine.disc_vars) == 3
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 1

    # Adaptive variable selection scheme :: disc_var_pick = 3
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 3,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "max_iter" => 2,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(JuMP.objective_value(m), 7049.2478976; atol = 1e-3)
    @test length(alpine.candidate_disc_vars) == 8
    @test length(alpine.disc_vars) == 8
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 3
end

@testset "Partitioning variable selection tests :: castro2m2" begin

    # Select all NL variables
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 0,
        "disc_uniform_rate" => 10,
        "presolve_bt" => false,
        "max_iter" => 0,
    )

    m = castro2m2(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED

    @test length(alpine.candidate_disc_vars) == 10
    @test length(alpine.disc_vars) == 10
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 0

    # Select minimum vertex cover
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 1,
        "disc_uniform_rate" => 10,
        "presolve_bt" => false,
        "max_iter" => 0,
    )

    m = castro2m2(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
    @test length(alpine.candidate_disc_vars) == 10
    @test length(alpine.disc_vars) == 4
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 1

    # Criteria 15 static selection
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 2,
        "disc_uniform_rate" => 15,
        "presolve_bt" => false,
        "max_iter" => 0,
    )

    m = castro2m2(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED

    @test length(alpine.candidate_disc_vars) == 10
    @test length(alpine.disc_vars) == 10
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 2
end

@testset "Partitioning variable selection tests :: blend029" begin

    # Select all NL variable
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => PAVITO,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 0,
        "disc_uniform_rate" => 10,
        "presolve_bt" => false,
        "max_iter" => 1,
        "log_level" => 100,
    )

    m = blend029_gl(solver = test_solver)
    alpine = _build(m)

    @test length(alpine.candidate_disc_vars) == 26
    @test Set(alpine.candidate_disc_vars) == Set([
        26,
        27,
        29,
        30,
        32,
        33,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        55,
        56,
        57,
        58,
        59,
        60,
    ])
    @test length(alpine.disc_vars) == 26
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 0

    # Minimum vertex cover
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => PAVITO,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 1,
        "disc_uniform_rate" => 10,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "max_iter" => 1,
        "log_level" => 100,
    )

    m = blend029_gl(solver = test_solver)
    alpine = _build(m)

    @test length(alpine.candidate_disc_vars) == 26
    @test Set(alpine.candidate_disc_vars) == Set([
        26,
        27,
        29,
        30,
        32,
        33,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        55,
        56,
        57,
        58,
        59,
        60,
    ])
    @test length(alpine.disc_vars) == 10
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 1

    # Adaptive Scheme vertex cover
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => PAVITO,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 2,
        "disc_uniform_rate" => 10,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "max_iter" => 1,
        "log_level" => 100,
    )

    m = blend029_gl(solver = test_solver)
    alpine = _build(m)

    @test length(alpine.candidate_disc_vars) == 26
    @test length(Set(alpine.candidate_disc_vars)) == 26
    # TODO provide a check to see if candidate_disc_vars are all covered
    @test length(alpine.disc_vars) == 10
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 2
end

@testset "Partitioning variable selection tests :: castro6m2" begin

    # Dynamic Scheme step 1
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_var_pick" => 3,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "max_iter" => 2,
        "log_level" => 100,
    )

    m = castro6m2(solver = test_solver)
    JuMP.optimize!(m)
    alpine = JuMP.backend(m).optimizer.model

    @test JuMP.termination_status(m) == MOI.OTHER_LIMIT
    @test JuMP.objective_value(m) <= 228.87

    @test length(alpine.candidate_disc_vars) == 24
    @test length(Set(alpine.candidate_disc_vars)) == 24
    @test length(alpine.disc_vars) == 12
    @test length(Set(alpine.disc_vars)) == 12
    @test MOI.get(m, MOI.RawOptimizerAttribute("disc_var_pick")) == 3
end

@testset "Hessians disabled with user-defined multivariate functions" begin
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "log_level" => 100,
    )
    m = Model(test_solver)
    my_f(x, y) = (x - 1)^2 + (y - 2)^2
    JuMP.register(m, :my_f, 2, my_f, autodiff = true)
    @variable(m, x[1:2])
    @NLobjective(m, Min, my_f(x[1], x[2]))

    JuMP.set_optimize_hook(m, MOI.Utilities.attach_optimizer)
    JuMP.optimize!(m)
    JuMP.set_optimize_hook(m, nothing)
    alpine = JuMP.backend(m).optimizer.model
    @test !(:Hess in Alpine.features_available(alpine))
end
