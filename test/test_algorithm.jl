@testset " Validation Test || AMP-TMC || basic solve || examples/nlp3.jl (2 iterations)" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => false,
        # "monomial_convexhull" => false,
        "presolve_bp" => true,
        "max_iter" => 2,
        "presolve_bt_width_tol" => 1e-3,
        "presolve_bt" => false,
        "partition_scaling_factor" => 10,
        "disc_var_pick" => 0,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT

    @test isapprox(objective_value(m), 7049.247897696512; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 2
end

@testset " Validation Test || AMP-TMC || minimum-vertex solving || examples/nlp3.jl (3 iterations)" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => false,
        # "monomial_convexhull" => false,
        "presolve_bp" => true,
        "disc_var_pick" => 1,
        "max_iter" => 2,
        "partition_scaling_factor" => 10,
        "presolve_bt_width_tol" => 1e-3,
        "presolve_bt" => false,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 7049.247897696512; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 2
end

@testset " Validation Test || BT-AMP-TMC || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => false,
        "log_level" => 100,
        "max_iter" => 2,
        "partition_scaling_factor" => 10,
        "presolve_bt_width_tol" => 1e-3,
        "presolve_bt_bound_tol" => 1e-1,
        "presolve_bt" => true,
        "presolve_bt_algo" => 1,
        "presolve_bp" => true,
        "presolve_bt_max_iter" => 2,
        "presolve_track_time" => true,
        "disc_var_pick" => max_cover_var_picker,
    )
    m = nlp3(solver = test_solver)

    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test MOI.get(m, Alpine.NumberOfIterations()) == 2
    @test MOI.get(m, Alpine.NumberOfPresolveIterations()) == 2

    vars = all_variables(m)
    LB = [100, 1000, 1000, 10, 150.2, 10, 35.4, 168]
    UB = [4573.6, 5547.8, 5913.3, 332.4, 551, 390, 571.1, 638.7]
    @test isapprox(MOI.get.(m, Alpine.TightenedLowerBound(), vars), LB, atol = 1E-6)
    @test isapprox(MOI.get.(m, Alpine.TightenedUpperBound(), vars), UB, atol = 1E-6)
end

# FIXME Pavito terminates with `NUMERICAL_ERROR` on Julia v1.0:
@testset " Validation Test || PBT-AMP-TMC || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => false,
        "log_level" => 100,
        "max_iter" => 2,
        "partition_scaling_factor" => 10,
        "presolve_bt" => true,
        "presolve_bt_width_tol" => 1e-3,
        "presolve_bt_bound_tol" => 1e-1,
        "presolve_bt_algo" => 2,
        "presolve_bp" => true,
        "presolve_bt_max_iter" => 2,
        "disc_var_pick" => max_cover_var_picker,
    )

    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    @test termination_status(m) == MOI.OTHER_LIMIT
    @test MOI.get(m, Alpine.NumberOfIterations()) == 2
end

@testset " Validation Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => true,
        "partition_scaling_factor" => 12,
        "max_iter" => 1,
    )
    m = nlp1(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 58.38367169858795; atol = 1e-5)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
end

@testset " Validation Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => false,
        "partition_scaling_factor" => 14,
        "max_iter" => 4,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 7049.24789; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 4
end

@testset " Validation Test || AMP || basic solve || examples/circle.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        "disc_abs_width_tol" => 1e-3,
        "partition_scaling_factor" => 8,
        "max_iter" => 6,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = circle(solver = test_solver)
    # This is fixed in Alpine - TODO (odow): cycling detected in Pavito when disc_abs_width_tol = 1E-2
    JuMP.optimize!(m)
    @test isapprox(objective_value(m), 1.4142135534556992; atol = 1e-3)
end

@testset " Validation Test || AMP-CONV-FACET || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => true,
        "convhull_formulation" => "facet",
        "max_iter" => 3,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test isapprox(objective_value(m), 7049.247897696188; atol = 1e-4)
    @test isapprox(objective_bound(m), 6654.6983279983; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 3
end

@testset " Validation Test || AMP || multi4N || N = 2 || exprmode=1:11" begin
    n_instances = 5
    objValVec = 2.0 * ones(n_instances)

    for i in 1:n_instances
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "disc_abs_width_tol" => 1e-2,
            "max_iter" => 4,
            "partition_scaling_factor" => 4,
            "presolve_bp" => false,
            "presolve_bt" => false,
            "log_level" => 1,
        )

        m = multi4N(solver = test_solver, N = 2, exprmode = i)
        JuMP.optimize!(m)

        @test termination_status(m) == MOI.OTHER_LIMIT
        @test isapprox(objective_value(m), objValVec[i]; atol = 1e-3)
    end
end

@testset " Validation Test || AMP || multi2 || exprmode=1:11" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "max_iter" => 5,
        "presolve_bt" => false,
        "log_level" => 1,
    )

    m = multi2(solver = test_solver)
    JuMP.optimize!(m)
    @test termination_status(m) == MOI.OPTIMAL
    @test primal_status(m) == MOI.FEASIBLE_POINT
    @test dual_status(m) == MOI.NO_SOLUTION
    @test raw_status(m) isa String
    @test isapprox(objective_value(m), 0.92906489; atol = 1e-3)
end

@testset " Validation Test || AMP || multi3N || N = 2 || exprmode=1:11" begin
    n_instances = 3
    objValVec = 2.0 * ones(n_instances)

    objBoundVec = Any[2.97186, 3.85492, 4.23375]

    for i in 1:n_instances
        test_solver = optimizer_with_attributes(
            Alpine.Optimizer,
            "nlp_solver" => IPOPT,
            "mip_solver" => HIGHS,
            "disc_abs_width_tol" => 1e-2,
            "max_iter" => 4,
            "partition_scaling_factor" => 4,
            "presolve_bp" => false,
            "presolve_bt" => false,
            "log_level" => 1,
        )

        m = multi3N(solver = test_solver, N = 2, exprmode = i)
        JuMP.optimize!(m)

        @test termination_status(m) == MOI.OTHER_LIMIT
        @test isapprox(objective_value(m), objValVec[i]; atol = 1e-3)
        @test objective_bound(m) <= objBoundVec[i] + 1e-3
    end
end

@testset " Validation Test || AMP || multiKND || K = 3, N = 3, D = 0 " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "max_iter" => 3,
        "partition_scaling_factor" => 4,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "log_level" => 1,
    )

    m = multiKND(solver = test_solver, randomub = 50, K = 3, N = 3, D = 0)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 3.0000000824779454; atol = 1e-3)
    @test isapprox(objective_bound(m), 12.054604248046875; atol = 1e-3)
end

@testset " Validation Test || AMP-CONV-FACET || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => false,
        "max_iter" => 3,
        "partition_scaling_factor" => 4,
        "convhull_formulation" => "facet",
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 7049.247897696188; atol = 1e-4)
    @test isapprox(objective_bound(m), 5871.530692199214; atol = 1E-5)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 3
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR || examples/nlp3.jl " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => false,
        "partition_scaling_factor" => 18,
        "max_iter" => 1,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 18
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/nlp3.jl " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 14
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/castro2m2.jl " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = castro2m2(solver = test_solver)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 8
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/multi3N.jl exprmode=2" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = multi3N(solver = test_solver, N = 3, exprmode = 1)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 16
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/multi3N.jl exprmode=2" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = multi3N(solver = test_solver, N = 3, exprmode = 1)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 20
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/multi4N.jl exprmode=1" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = multi4N(solver = test_solver, N = 2, exprmode = 1)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 12
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/multi4N.jl exprmode=2" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => false,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = multi4N(solver = test_solver, N = 2, exprmode = 1)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 20
end

@testset " Validation Test || AMP || PARTITION-SCALING-FACTOR-BRANCH || examples/multi4N.jl exprmode=2" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_abs_width_tol" => 1e-2,
        "partition_scaling_factor_branch" => true,
        "max_iter" => 1,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = multi4N(solver = test_solver, N = 2, exprmode = 2)
    JuMP.optimize!(m)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test MOI.get(m, MOI.RawOptimizerAttribute("partition_scaling_factor")) == 20
end

@testset "Operator :: bmpl && binlin && binprod solve test I" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => PAVITO,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => false,
        "log_level" => 100,
    )

    m = bpml_lnl(solver = test_solver)
    JuMP.optimize!(m)
    @test isapprox(objective_value(m), 0.3; atol = 1e-6)
    alpine = JuMP.backend(m).optimizer.model
    @test haskey(alpine.nonconvex_terms, Expr[:(x[1]), :(x[6])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[2]), :(x[7])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[3]), :(x[8])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[4]), :(x[9])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[5]), :(x[10])])

    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:nonlinear_type] == :BINLIN
end

@testset "Operator :: bmpl && binlin && binprod solve test II" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => JUNIPER,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => false,
        "partition_scaling_factor" => 4,
        "log_level" => 100,
    )

    m = bpml_binl(solver = test_solver)

    # FIXME: Deactivating this until Juniper v0.9.0's numerical issues are fixed. 
    # JuMP.optimize!(m)
    # @test isapprox(JuMP.objective_value(m), 22812.76415926; atol=1e-6)
    # alpine = JuMP.backend(m).optimizer.model

    alpine = _build(m)
    @test haskey(alpine.nonconvex_terms, Expr[:(x[6]), :(x[7])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[7]), :(x[8])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[8]), :(x[9])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[9]), :(x[10])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[10]), :(x[6])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[1]), :(x[11])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[2]), :(x[13])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[3]), :(x[15])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[4]), :(x[17])])
    @test haskey(alpine.nonconvex_terms, Expr[:(x[5]), :(x[19])])

    @test alpine.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[13])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[15])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[4]), :(x[17])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[5]), :(x[19])]][:nonlinear_type] == :BINLIN
end

@testset "Embedding Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => true,
        "convhull_ebd" => true,
        "max_iter" => 1,
        "partition_scaling_factor" => 4,
        "log_level" => 100,
    )

    m = nlp1(solver = test_solver)
    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 58.38367169858795; atol = 1e-4)
    @test isapprox(alp.best_bound, 52.9702, atol = 1E-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
end

@testset "Embedding IBS Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => true,
        "convhull_ebd" => true,
        "convhull_ebd_ibs" => true,
        "partition_scaling_factor" => 8,
        "max_iter" => 1,
        "log_level" => 100,
    )
    m = nlp1(solver = test_solver)
    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 58.38367169858795; atol = 1e-5)
    @test isapprox(alp.best_bound, 57.012, atol = 1E-4)

    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
end

@testset "Embedding IBS Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => false,
        "convhull_ebd" => true,
        "convhull_ebd_ibs" => true,
        "max_iter" => 3,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 7049.247897696188; atol = 1e-5)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 3
    @test isapprox(alp.best_bound, 6893.2067, atol = 1E-4)
end

@testset "Embedding IBS Test || AMP || special problem || ... " begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        "disc_abs_width_tol" => 1e-3,
        "partition_scaling_factor" => 8,
        "max_iter" => 6,
        "presolve_bt" => false,
        "presolve_bp" => true,
        "presolve_bt_algo" => 1,
        "presolve_bt_bound_tol" => 1e-1,
        "convhull_ebd" => true,
        "convhull_ebd_ibs" => true,
        "log_level" => 100,
    )

    m = circle(solver = test_solver)
    # This is fixed in Alpine - TODO(odow): mixed-integer cycling detected, terminating Pavito
    JuMP.optimize!(m)
    @test isapprox(objective_value(m), 1.41421355; atol = 1e-3)
end

@testset "Embedding LINK Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bt" => false,
        "presolve_bp" => false,
        "convhull_ebd" => true,
        "convhull_ebd_link" => true,
        "partition_scaling_factor" => 10,
        "max_iter" => 3,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 7049.247897696188; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 3
end

@testset "Algorithm Logic Test || castro4m2 || 1 iteration || Error case" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "max_iter" => 1,
        "presolve_bt" => false,
        "partition_scaling_factor" => 10,
        "log_level" => 100,
    )

    m = castro4m2(solver = test_solver)
    JuMP.optimize!(m)
    @test termination_status(m) == MOI.OTHER_LIMIT
end

@testset " Algorithm Logic Test || blend029_gl || 3 iterations || Infeasible Case" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => PAVITO,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bp" => true,
        "disc_var_pick" => 1,
        "log_level" => 100,
        "max_iter" => 3,
        "partition_scaling_factor" => 10,
        "presolve_bt_width_tol" => 1e-3,
        "presolve_bt" => false,
    )
    m = blend029_gl(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test MOI.get(m, Alpine.NumberOfIterations()) == 3
    @test objective_bound(m) <= 14.0074
end

@testset "Convex Model Solve" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        "max_iter" => 1,
        "presolve_bt" => false,
        "log_level" => 100,
    )
    m = convex_solve(solver = test_solver)
    JuMP.optimize!(m)
    @test termination_status(m) == MOI.OPTIMAL
end

@testset "Uniform partitioning" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "disc_add_partition_method" => "uniform",
        "disc_uniform_rate" => 10,
        "max_iter" => 1,
        "presolve_bt" => false,
        "time_limit" => 100000,
        "log_level" => 100,
    )
    m = nlp3(solver = test_solver)
    JuMP.optimize!(m)
    @test isapprox(objective_bound(m), 6561.7841; atol = 1e-3)
end

@testset "Algorithm Test with binprod terms" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => PAVITO,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        # "bilinear_convexhull" => true,
        # "monomial_convexhull" => true,
        "presolve_bp" => true,
        "presolve_bt" => false,
        "partition_scaling_factor" => 10,
        "log_level" => 100,
    )
    m = binprod_nlp3(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OPTIMAL
    @test isapprox(objective_value(m), 3651.020370626844; rtol = 1e-4)
    @test isapprox(objective_bound(m), 3650.786358471944; rtol = 1e-4)

    alpine = JuMP.backend(m).optimizer.model
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:y_idx] == 19
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:id] == 6
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:lifted_constr_ref] ==
          :(x[19] == x[2] * x[4])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:y_idx] == 25
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:id] == 12
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:lifted_constr_ref] ==
          :(x[25] == x[3] * x[8])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:y_idx] == 22
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:id] == 9
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:lifted_constr_ref] ==
          :(x[22] == x[3] * x[5])
    @test alpine.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:y_idx] == 20
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:id] == 7
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:lifted_constr_ref] ==
          :(x[20] == x[12] * x[19])
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:y_idx] == 14
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:id] == 1
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:lifted_constr_ref] ==
          :(x[14] == x[9] * x[4])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:y_idx] == 17
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:id] == 4
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:lifted_constr_ref] ==
          :(x[17] == x[1] * x[6])
    @test alpine.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:y_idx] == 16
    @test alpine.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:id] == 3
    @test alpine.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:lifted_constr_ref] ==
          :(x[16] == x[11] * x[5])
    @test alpine.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:y_idx] == 30
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:id] == 17
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:lifted_constr_ref] ==
          :(x[30] == x[10] * x[11])
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:nonlinear_type] == :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:y_idx] == 29
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:id] == 16
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:lifted_constr_ref] ==
          :(x[29] == x[9] * x[10] * x[11])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:nonlinear_type] ==
          :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:y_idx] == 21
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:id] == 8
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:lifted_constr_ref] ==
          :(x[21] == x[2] * x[7])
    @test alpine.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:y_idx] == 32
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:id] == 19
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:lifted_constr_ref] ==
          :(x[32] == x[9] * x[13])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:nonlinear_type] == :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:y_idx] == 15
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:id] == 2
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:lifted_constr_ref] ==
          :(x[15] == x[10] * x[6])
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:y_idx] == 33
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:id] == 20
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:lifted_constr_ref] ==
          :(x[33] == x[10] * x[12])
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:nonlinear_type] == :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:y_idx] == 28
    @test alpine.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:id] == 15
    @test alpine.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:lifted_constr_ref] ==
          :(x[28] == x[27] * x[5])
    @test alpine.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:y_idx] == 26
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:id] == 13
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:lifted_constr_ref] ==
          :(x[26] == x[13] * x[25])
    @test alpine.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:y_idx] == 27
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:id] == 14
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:lifted_constr_ref] ==
          :(x[27] == x[9] * x[12])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:nonlinear_type] == :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:y_idx] == 23
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:id] == 10
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:lifted_constr_ref] ==
          :(x[23] == x[10] * x[13])
    @test alpine.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:nonlinear_type] == :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:y_idx] == 24
    @test alpine.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:id] == 11
    @test alpine.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:lifted_constr_ref] ==
          :(x[24] == x[23] * x[22])
    @test alpine.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:nonlinear_type] == :BINLIN
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:y_idx] == 31
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:id] == 18
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:lifted_constr_ref] ==
          :(x[31] == x[12] * x[13])
    @test alpine.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:nonlinear_type] == :BINPROD
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:y_idx] == 18
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:id] == 5
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:lifted_constr_ref] ==
          :(x[18] == x[9] * x[17])
    @test alpine.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:nonlinear_type] == :BINLIN
    @test alpine.bounding_constr_mip[1][:rhs] == 1.0
    @test alpine.bounding_constr_mip[1][:vars] == Any[:(x[14]), :(x[15])]
    @test alpine.bounding_constr_mip[1][:coefs] == Any[0.0025, 0.0025]
    @test alpine.bounding_constr_mip[1][:sense] == :(<=)
    @test alpine.bounding_constr_mip[1][:cnt] == 2
    @test alpine.bounding_constr_mip[2][:rhs] == 1.0
    @test alpine.bounding_constr_mip[2][:vars] == Any[:(x[14]), :(x[5]), :(x[7])]
    @test alpine.bounding_constr_mip[2][:coefs] == Any[-0.0025, 0.0025, 0.0025]
    @test alpine.bounding_constr_mip[2][:sense] == :(<=)
    @test alpine.bounding_constr_mip[2][:cnt] == 3
    @test alpine.bounding_constr_mip[3][:rhs] == 1.0
    @test alpine.bounding_constr_mip[3][:vars] == Any[:(x[16]), :(x[8])]
    @test alpine.bounding_constr_mip[3][:coefs] == Any[-0.01, 0.01]
    @test alpine.bounding_constr_mip[3][:sense] == :(<=)
    @test alpine.bounding_constr_mip[3][:cnt] == 2
    @test alpine.bounding_constr_mip[4][:rhs] == 83333.333
    @test alpine.bounding_constr_mip[4][:vars] == Any[:(x[1]), :(x[18]), :(x[14])]
    @test alpine.bounding_constr_mip[4][:sense] == :(<=)
    @test alpine.bounding_constr_mip[4][:cnt] == 3
    @test alpine.bounding_constr_mip[5][:rhs] == 0.0
    @test alpine.bounding_constr_mip[5][:vars] ==
          Any[:(x[20]), :(x[21]), :(x[4]), :(x[5])]
    @test alpine.bounding_constr_mip[5][:coefs] == Any[1.0, -1.0, -1250.0, 1250.0]
    @test alpine.bounding_constr_mip[5][:sense] == :(<=)
    @test alpine.bounding_constr_mip[5][:cnt] == 4
    @test alpine.bounding_constr_mip[6][:rhs] == -1.25e6
    @test alpine.bounding_constr_mip[6][:vars] == Any[:(x[24]), :(x[26]), :(x[28])]
    @test alpine.bounding_constr_mip[6][:coefs] == Any[1.0, -1.0, -2500.0]
    @test alpine.bounding_constr_mip[6][:sense] == :(<=)
    @test alpine.bounding_constr_mip[6][:cnt] == 3
    @test alpine.bounding_constr_mip[7][:rhs] == 0.0
    @test alpine.bounding_constr_mip[7][:vars] == Any[:(x[29])]
    @test alpine.bounding_constr_mip[7][:coefs] == Any[1.0]
    @test alpine.bounding_constr_mip[7][:sense] == :(<=)
    @test alpine.bounding_constr_mip[7][:cnt] == 1
    @test alpine.bounding_constr_mip[8][:rhs] == 0.0
    @test alpine.bounding_constr_mip[8][:vars] == Any[:(x[30]), :(x[31])]
    @test alpine.bounding_constr_mip[8][:coefs] == Any[1.0, -1.0]
    @test alpine.bounding_constr_mip[8][:sense] == :(>=)
    @test alpine.bounding_constr_mip[8][:cnt] == 2
    @test alpine.bounding_constr_mip[9][:rhs] == 0.0
    @test alpine.bounding_constr_mip[9][:vars] == Any[:(x[32]), :(x[33])]
    @test alpine.bounding_constr_mip[9][:coefs] == Any[1.0, -1.0]
    @test alpine.bounding_constr_mip[9][:sense] == :(<=)
    @test alpine.bounding_constr_mip[9][:cnt] == 2
end

@testset "TESTS for closing the optimality gap in OBBT" begin
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => PAVITO,
        "presolve_bt" => true,
        "presolve_bt_max_iter" => 2,
        "log_level" => 1,
    )
    m = JuMP.Model(test_solver)

    # From issue #108
    @variable(m, -2 <= x[1:2] <= 2)
    @variable(m, -10 <= y[1:3] <= 10)
    @NLobjective(m, Min, y[2] + y[1] * y[2] * y[3])
    @constraint(m, y[1] == x[1])
    @constraint(m, y[2] == x[2])
    @NLconstraint(m, y[3] == x[2]^2)

    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model
    @test isapprox(JuMP.objective_value(m), -18, atol = 1E-6)
    @test isapprox(value.(m[:x]), [2, -2], atol = 1E-6)
    @test alp.logs[:n_iter] == 0
    @test MOI.get(m, Alpine.NumberOfIterations()) == 0
end

@testset "Linking constraints for multilinear terms" begin

    # Without linking constraints
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => false,
        "apply_partitioning" => true,
        "linking_constraints" => false,
    )
    m = linking_constraints_testing(solver = test_solver)
    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model
    @test isapprox(JuMP.objective_value(m), -0.5294000135, atol = 1E-5)
    @test isnothing(alp.linking_constraints_info)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 2

    # With linking constraints
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => false,
        "apply_partitioning" => true,
        "linking_constraints" => true,
    )
    m = linking_constraints_testing(solver = test_solver)
    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model
    @test isapprox(JuMP.objective_value(m), -0.5294000135, atol = 1E-5)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 1
    @test in(Set{Any}([1, 2]), alp.linking_constraints_info[[1, 2]])
    @test in(Set{Any}([1, 2, 3]), alp.linking_constraints_info[[1, 2]])
    @test in(Set{Any}([1, 2, 4]), alp.linking_constraints_info[[1, 2]])
    @test in(Set{Any}([1, 3]), alp.linking_constraints_info[[1, 3]])
    @test in(Set{Any}([1, 2, 3]), alp.linking_constraints_info[[1, 3]])
    @test in(Set{Any}([1, 3, 4]), alp.linking_constraints_info[[1, 3]])
    @test in(Set{Any}([1, 4]), alp.linking_constraints_info[[1, 4]])
    @test in(Set{Any}([1, 2, 4]), alp.linking_constraints_info[[1, 4]])
    @test in(Set{Any}([1, 3, 4]), alp.linking_constraints_info[[1, 4]])
end

@testset "Warm start value used as the incumbent solution" begin
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => false,
        "max_iter" => 1,
        "use_start_as_incumbent" => false,
    )
    sol = [
        579.30669,
        1359.97066,
        5109.97054,
        182.0177,
        295.60118,
        217.9823,
        286.41653,
        395.60118,
    ]

    m = nlp3(solver = test_solver)
    for i in 1:length(m[:x])
        JuMP.set_start_value(m[:x][i], sol[i])
    end

    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(objective_value(m), 7049.247897696512, atol = 1e-6)
    @test isapprox(alp.best_bound, 3834.978106740535, atol = 1E-6)
    @test alp.options.use_start_as_incumbent == false

    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "presolve_bt" => true,
        "max_iter" => 1,
        "use_start_as_incumbent" => true,
    )

    m = nlp3(solver = test_solver)
    for i in 1:length(m[:x])
        JuMP.set_start_value(m[:x][i], sol[i])
    end

    JuMP.optimize!(m)
    alp = JuMP.backend(m).optimizer.model

    @test termination_status(m) == MOI.OTHER_LIMIT
    @test isapprox(alp.best_bound, 4810.212866711817, atol = 1E-3)
    @test alp.warm_start_value == sol
    @test alp.options.use_start_as_incumbent == true
end

@testset "Test binary multilinear product linearization" begin
    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "minlp_solver" => JUNIPER,
        "presolve_bt" => false,
    )

    m = hmittelman(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OPTIMAL
    @test isapprox(objective_value(m), 13; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 1

    test_solver = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "minlp_solver" => JUNIPER,
        "presolve_bt" => true,
    )

    m = hmittelman(solver = test_solver)
    JuMP.optimize!(m)

    @test termination_status(m) == MOI.OPTIMAL
    @test isapprox(objective_value(m), 13; atol = 1e-4)
    @test MOI.get(m, Alpine.NumberOfIterations()) == 0
end

@testset "Test integer variable support via IntegerToZeroOneBridge" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "minlp_solver" => JUNIPER,
    )
    model = JuMP.Model(test_solver)
    @variable(model, -10 <= x <= 20, Int)
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test isapprox(objective_value(model), -10; atol = 1e-4)
    @test isapprox(value(x), -10, atol = 1e-6)
    @test MOI.get(model, Alpine.NumberOfIterations()) == 0
end

@testset "Test integer variable support 0" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "minlp_solver" => JUNIPER,
    )
    model = JuMP.Model(test_solver)
    @variable(model, x, Int)
    @objective(model, Min, x)
    @test_throws(
        ErrorException(
            "Unable to use IntegerToZeroOneBridge because the variable MOI.VariableIndex(1) has a non-finite domain",
        ),
        JuMP.optimize!(model),
    )
end

@testset "Test integer variable support 1" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "minlp_solver" => JUNIPER,
    )
    model = JuMP.Model(test_solver)
    @variable(model, -10 <= x, Int)
    @objective(model, Min, x)
    @test_throws(
        ErrorException(
            "Unable to use IntegerToZeroOneBridge because the variable MOI.VariableIndex(1) has a non-finite domain",
        ),
        JuMP.optimize!(model),
    )
end

@testset "Test integer variable support 2" begin
    test_solver = optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => IPOPT,
        "mip_solver" => HIGHS,
        "minlp_solver" => JUNIPER,
    )
    model = JuMP.Model(test_solver)
    @variable(model, x <= 20, Int)
    @objective(model, Min, x)
    @test_throws(
        ErrorException(
            "Unable to use IntegerToZeroOneBridge because the variable MOI.VariableIndex(1) has a non-finite domain",
        ),
        JuMP.optimize!(model),
    )
end
