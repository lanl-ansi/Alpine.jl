
@testset "linear constraints FBBT tests" begin
    solver_options = DefaultTestSolver(log_level=0)
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

    optimize!(m)

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


end 