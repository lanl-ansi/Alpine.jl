
@testset "loading tests" begin
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
end 