
@testset "loading/printing tests" begin
    solver_options = DefaultTestSolver(log_level=0)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )

    @variable(m, -5 <= x[1:2] <= 5)
    @objective(m, Min, 6*x[1]^2 + 4*x[2]^2 - 2.5*x[1]*x[2])
    @NLconstraint(m, x[1]*x[2] >= 8)
    optimize!(m)

    bm = JuMP.backend(m)
    internal_model = bm.optimizer.model

    @test internal_model.inner.is_objective_quadratic == true 
    @test length(internal_model.inner.nl_constraint_expr) == 1
    @test length(internal_model.variable_info) == 2
    @test internal_model.inner.objective_convexity == :convex

    m = Model(with_optimizer(
        Alpine.Optimizer, 
        nlp_optimizer = nlp_optimizer(), 
        mip_optimizer = mip_optimizer()
    ))

    @variable(m, x, Bin)
    @variable(m, y, start = 0.0)
    JuMP.set_upper_bound(y, 5)
    JuMP.set_lower_bound(y, -5)
    @constraint(m, y == 0.0)
    @objective(m, Min, y^2)
    @constraint(m, x * y <= 3)
    optimize!(m)

    bm = JuMP.backend(m)
    internal_model = bm.optimizer.model

    @test internal_model.inner.num_constraints == 2 
    @test internal_model.inner.is_objective_quadratic == true
    @test internal_model.inner.objective_convexity == :convex 
    @test internal_model.inner.quadratic_function_convexity[:quadratic_le][1] == :undet
    @test internal_model.inner.quadratic_constraint_convexity[:quadratic_le][1] == :undet
 
end 