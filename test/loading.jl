
@testset "model loading tests" begin
    
    printstyled("\nTesting model loading ...\n", color=:blue, bold=true) 

    solver_options = DefaultTestSolver(log_level=0, 
        perform_bp_only=true, 
        max_multistart_points=0)
    
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )

    @variable(m, -5 <= x[1:2] <= 5)
    @objective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5*x[1]*x[2])
    @NLconstraint(m, x[1]*x[2] >= 8)
    @constraint(m, x[1]^2 <= -x[2]^2 + 100)
    @constraint(m, x[1]^2 == 25)
    JuMP.optimize!(m)

    bm = JuMP.backend(m)
    internal_model = bm.optimizer.model

    @test internal_model.inner.is_objective_quadratic == true 
    @test length(internal_model.inner.nl_constraint_expression) == 1
    @test length(internal_model.variable_info) == 2

    solver_options = DefaultTestSolver(log_level=0, 
        perform_bp_only=true, 
        max_multistart_points=0)
    
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )

    @variable(m, x, Bin)
    @variable(m, y, start = 0.0)
    JuMP.set_upper_bound(y, 5)
    JuMP.set_lower_bound(y, -5)
    @constraint(m, y == 0.0)
    @objective(m, Min, y^2)
    @constraint(m, x * y <= 3)
    JuMP.optimize!(m)

    bm = JuMP.backend(m)
    internal_model = bm.optimizer.model

    @test internal_model.inner.num_constraints == 2 
    @test internal_model.inner.is_objective_quadratic == true

    # trig instance from MINLPLib2
    solver_options = DefaultTestSolver(log_level=0, 
        perform_bp_only=true, 
        max_multistart_points=100)
    
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )

    @variable(m, -2 <= x <= 5)
    JuMP.set_start_value(x, 1.0)
    
    @NLconstraint(m, 5*sin(x) - x <= 0)
    @NLobjective(m, Min, sin(11*x) + cos(13*x) - sin(17*x) - cos(19*x)) 
    JuMP.optimize!(m)

    @test isapprox(JuMP.objective_value(m), -3.7625, atol=1e-1)
    
end 