
@testset "FBBT tests" begin

    printstyled("\nTesting FBBT ...\n", color=:blue, bold=true)

    solver_options = DefaultTestSolver(log_level=0, perform_bp_only=true)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )

    @variable(m, x[1:2])
    JuMP.set_lower_bound(x[1], 0)
    JuMP.set_upper_bound(x[1], 1)
    @constraint(m, sum(x) == 0)

    @variable(m, y[1:2])
    JuMP.set_lower_bound(y[1], 0)
    JuMP.set_upper_bound(y[1], 1)
    @constraint(m, sum(y) <= 0)
    @constraint(m, sum(y) >= 0)

    JuMP.optimize!(m)

    bm = JuMP.backend(m)
    internal_model = bm.optimizer.model
    alpine_problem = internal_model.inner

    @test alpine_problem.variable_bound_original[2].hi ==  Inf 
    @test alpine_problem.variable_bound_original[2].lo == -Inf
    @test alpine_problem.variable_bound_tightened[2].hi == 0
    @test alpine_problem.variable_bound_tightened[2].lo == -1

    @test alpine_problem.variable_bound_original[4].hi ==  Inf 
    @test alpine_problem.variable_bound_original[4].lo == -Inf
    @test alpine_problem.variable_bound_tightened[4].hi == 0
    @test alpine_problem.variable_bound_tightened[4].lo == -1

    # infeasible LP 
    solver_options = DefaultTestSolver(log_level=0, perform_bp_only=true)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )
    @variable(m, 0 <= x[1:2] <= 1)
    @constraint(m, sum(x) >= 3)

    JuMP.optimize!(m) 

    @test JuMP.termination_status(m) == MOI.INFEASIBLE

    # infeasible QP 
    solver_options = DefaultTestSolver(log_level=0, perform_bp_only=true)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )
    @variable(m, x[1:2])
    @constraint(m, x[1]^2 + x[2]^2 <= -1)

    JuMP.optimize!(m) 

    @test JuMP.termination_status(m) == MOI.INFEASIBLE

    # tighten bounds 
    solver_options = DefaultTestSolver(log_level=0, perform_bp_only=true)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options
    ))
    @variable(m, x)
    @variable(m, y)
    @variable(m, 0 <= z <= 1)
    @NLconstraint(m, exp(x^2 + y^2)-z <= 0)
    @NLconstraint(m, z*sqrt(x^2 + y^2) <= 1)

    JuMP.optimize!(m)

    lb_x = Alpine.lower_bound_tightened(x)
    ub_x = Alpine.upper_bound_tightened(x)

    @test isthin(lb_x..ub_x)
    @test isthin(Alpine.variable_bound_tightened(y))
    @test isthin(Alpine.variable_bound_tightened(z))

    # nvs01 from MINLPLib with continuous relaxation with changed bounds
    solver_options = DefaultTestSolver(log_level=0, perform_bp_only=true)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options 
    ))
    
    @variable(m, x[1:3])
    @variable(m, x_obj)

    JuMP.set_lower_bound(x[1], 200)
    JuMP.set_lower_bound(x[2], 200)
    JuMP.set_lower_bound(x[3], 0)
    JuMP.set_upper_bound(x[1], 200)
    JuMP.set_upper_bound(x[2], 200)
    JuMP.set_upper_bound(x[3], 100)

    @NLconstraint(m, 420.169404664517*sqrt(900 + x[2]^2) - x[3]*x[1]*x[2] == 0)
    @NLconstraint(m, (2960.87631843 + 18505.4769901875*x[2]^2)/(7200 + x[1]^2) - x[3] >= 0)
    @NLconstraint(m, -0.04712385*sqrt(900 + x[1]^2)*x[2] + x_obj >= 0)

    JuMP.optimize!(m)

    @test isapprox(Alpine.lower_bound_tightened(x_obj), 1906.04, atol=1e-1)
    @test isapprox(Alpine.lower_bound_tightened(x[3]), 2.12435, atol=1e-1)
    @test isapprox(Alpine.upper_bound_tightened(x[3]), 2.12435, atol=1e-1)

    # redundant constraints 
    solver_options = DefaultTestSolver(perform_bp_only=true)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options 
    ))

    @variable(m, 0 <= x[1:3] <= 1)

    @constraint(m, sum(x) <= 100)
    @constraint(m, x[1]^2 + x[2]^2 + x[3]^2 >= -10)
    @constraint(m, x[1]^2 + x[2]^2 + x[3]^2 <= 10)
    @NLconstraint(m, x[1]^3 + x[2]^3 <= 10)

    JuMP.optimize!(m)
    bm = JuMP.backend(m)
    internal_model = bm.optimizer.model
    alpine_problem = internal_model.inner

    @test length(alpine_problem.redundant_constraint_ids.linear_constraint_ids) == 1
    @test length(alpine_problem.redundant_constraint_ids.quadratic_constraint_ids) == 2
    @test length(alpine_problem.redundant_constraint_ids.nl_constraint_ids) == 1

end 