
include("instances/gtm.jl")
@testset "convexity detection tests" begin

    printstyled("\nTesting quadratic function convexity detection ...\n", color=:blue, bold=true)    
    Alpine.silence()
    optimizer = nlp_optimizer()
    num_variables = 4
    x = MOI.add_variables(optimizer, num_variables)
    
    quadratic_terms = []
    affine_terms = Alpine.SAF[]
    constant = 0.0
    push!(quadratic_terms, Alpine.SQT(1.0, x[1], x[1]))
    push!(quadratic_terms, Alpine.SQT(1.0, x[2], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :convex

    push!(quadratic_terms, Alpine.SQT(-2.0, x[1], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :convex
    
    quadratic_terms = []
    push!(quadratic_terms, Alpine.SQT(-1.0, x[1], x[1]))
    push!(quadratic_terms, Alpine.SQT(-1.0, x[2], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :concave

    quadratic_terms = []
    push!(quadratic_terms, Alpine.SQT(-1.0, x[1], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :undet

    # convex constraint detection 
    solver_options = DefaultTestSolver(log_level=1)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options 
    ))

    @variable(m, 2 <= x <= 3)
    @variable(m, 2 <= y <= 3)
    @variable(m, 0 <= z <= 1)
    
    @NLconstraint(m, (log(x+y-3*z^2))^(-2) <= 10)

    JuMP.optimize!(m)


    @test Alpine.is_problem_convex(m) == true
    @test isapprox(JuMP.objective_value(m), 0.0, atol=1e-1)
    
    # convex problem detection
    solver_options = DefaultTestSolver(log_level=1)
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options 
    ))
    
    m = gtm(m = m)

    JuMP.optimize!(m)

    @test Alpine.is_problem_convex(m) == true
    @test isapprox(JuMP.objective_value(m), 543.5651, atol=1e-1)

    
end 