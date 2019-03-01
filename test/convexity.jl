
@testset "function convexity detection" begin
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

    
end 