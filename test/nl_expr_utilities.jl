
@testset "nonlinear expression utilities" begin
    
    solver_options = DefaultTestSolver(log_level=0)
    @test MOI.supports(Alpine.Optimizer(solver_options), MOI.NLPBlock()) == true
    
    Alpine.silence()
    optimizer = nlp_optimizer()
    num_variables = 3
    x = MOI.add_variables(optimizer, num_variables)

    expressions = Expr[
        :(x[$(x[1])] + (0.5 + 0.25 * 1.0)        == 2.0), # constant sub-tree flattening
        :(x[$(x[1])] * -x[$(x[2])] * +x[$(x[3])] == 2.0), # multilinear sign separation
        :(-2.0 * -x[$(x[2])] * -4.0              == 2.0), # multilinear coeff aggregation
    ]

    Alpine.expr_flatten_constant_subtree(expressions[1])
    Alpine.expr_separate_sign_multilinear(expressions[2])
    Alpine.expr_separate_sign_multilinear(expressions[3])
    Alpine.expr_aggregate_coeff_multilinear(expressions[3])
    
    @test expressions[1] == :(x[$(x[1])] + 0.75 == 2.0)
    @test expressions[2] == :(x[$(x[1])] * x[$(x[2])] * x[$(x[3])] * -1.0 == 2.0)
    @test expressions[3] == :(-8.0 * x[$(x[2])] == 2.0)
end