
@testset "nonlinear expression utilities" begin
    
    solver_options = DefaultTestSolver(log_level=0)
    @test MOI.supports(Alpine.Optimizer(solver_options), MOI.NLPBlock()) == true
    
    Alpine.silence()
    optimizer = nlp_optimizer()
    num_variables = 3
    x = MOI.add_variables(optimizer, num_variables)

    expressions = Expr[
        :(x[$(x[1])] + (0.5 + 0.25 * 1.0) == 2.0), # 1
        :(x[$(x[1])] * -x[$(x[2])] * +x[$(x[3])] == 2.0), # 2
        :(-2.0 * -x[$(x[2])] * -4.0 == 2.0), # 3
        :(x[$(x[1])] - (sin(x[$(x[2])]) - -4.0 * log(x[$(x[3])]) *x[$(x[1])]) == 2.0), # 4
        :(x[$(x[1])] + x[$(x[1])] - x[$(x[1])]^2 + x[$(x[1])] * x[$(x[2])] - 
        x[$(x[1])]^2 + x[$(x[2])] * x[$(x[1])] + x[$(x[2])] * x[$(x[1])] * x[$(x[3])] +
        x[$(x[1])] * x[$(x[3])] * x[$(x[2])] == 2.0), # 5
        :(x[$(x[1])]^5 + 2.5 * x[$(x[1])]^5 == 2.0),  # 6
        :(log(x[$(x[1])]) + 5.0 * log(x[$(x[1])]) == 2.0), # 7
        :(sin(x[$(x[1])]) + cos(x[$(x[2])]) + 3.0 * sin(x[$(x[1])]) == 2.0), # 8
        :(abs(x[$(x[1])]) + 2.0 * abs(x[$(x[1])]) == 2.0), # 9
        :(log(x[$(x[1])]) * exp(x[$(x[1])]) == 2.0), # 10
        :(exp(x[$(x[1])]) + 5.0 * exp(x[$(x[1])]) == 2.0), # 11
    ]

    Alpine.expr_flatten_constant_subtree(expressions[1])
    Alpine.expr_separate_sign_multilinear(expressions[2])
    Alpine.expr_separate_sign_multilinear(expressions[3])
    Alpine.expr_aggregate_coeff_multilinear(expressions[3])
    
    @test expressions[1] == :(x[$(x[1])] + 0.75 == 2.0)
    @test expressions[2] == :(x[$(x[1])] * x[$(x[2])] * x[$(x[3])] * -1.0 == 2.0)
    @test expressions[3] == :(-8.0 * x[$(x[2])] == 2.0)

    disaggretated_expr = Alpine.expr_disaggregate(expressions[4])
    @test disaggretated_expr[1] ==  ( 1.0, :(x[$(x[1])]))  
    @test disaggretated_expr[2] ==  (-1.0, :(sin(x[$(x[2])])))
    @test disaggretated_expr[3] ==  (-4.0, :(log(x[$(x[3])]) * x[$(x[1])]))

    disaggretated_expr = Alpine.expr_disaggregate(expressions[5])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 8 
    
    @test nl_function.linear_part[1].expression[1] == 2.0 
    @test nl_function.linear_part[1].expression[2] == :(x[$(x[1])])
    @test nl_function.linear_part[1].convexity == :convex
    
    @test nl_function.quadratic_part[1].expression[1] == -2.0
    @test nl_function.quadratic_part[1].expression[2] == :(x[$(x[1])]^2)
    @test nl_function.quadratic_part[1].convexity == :concave

    @test nl_function.bilinear_part[1].expression[1] == 2.0
    @test nl_function.bilinear_part[1].expression[2] == :(x[$(x[1])] * x[$(x[2])])
    @test nl_function.bilinear_part[1].convexity == :undet
    
    @test nl_function.multilinear_part[1].expression[1] == 2.0
    @test nl_function.multilinear_part[1].expression[2] == :(x[$(x[1])] * x[$(x[2])] * x[$(x[3])])
    @test nl_function.multilinear_part[1].convexity == :undet

    disaggretated_expr = Alpine.expr_disaggregate(expressions[6])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 2

    @test nl_function.power_part[1].expression[1] == 3.5 
    @test nl_function.power_part[1].expression[2] == :(x[$(x[1])]^5)
    @test nl_function.power_part[1].convexity == :undet 

    disaggretated_expr = Alpine.expr_disaggregate(expressions[7])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 2 

    @test nl_function.log_part[1].expression[1] == 6.0
    @test nl_function.log_part[1].expression[2] == :(log(x[$(x[1])]))
    @test nl_function.log_part[1].convexity == :undet 

    disaggretated_expr = Alpine.expr_disaggregate(expressions[8])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 3 

    @test nl_function.trigonometric_part[1].expression[1] == 4.0
    @test nl_function.trigonometric_part[2].expression[1] == 1.0

    disaggretated_expr = Alpine.expr_disaggregate(expressions[9])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 2 

    @test nl_function.abs_part[1].expression[1] == 3.0
    @test nl_function.abs_part[1].expression[2] == :(abs(x[$(x[1])]))
    @test nl_function.abs_part[1].convexity == :undet
    
    disaggretated_expr = Alpine.expr_disaggregate(expressions[10])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 1 

    @test nl_function.other_part[1].expression[1] == 1.0
    @test nl_function.other_part[1].expression[2] == :(log(x[$(x[1])]) * exp(x[$(x[1])]))
    @test nl_function.other_part[1].convexity == :undet 

    disaggretated_expr = Alpine.expr_disaggregate(expressions[11])
    nl_function = Alpine.create_nl_function(disaggretated_expr)

    @test length(disaggretated_expr) == 2 

    @test nl_function.exp_part[1].expression[1] == 6.0
    @test nl_function.exp_part[1].expression[2] == :(exp(x[$(x[1])]))
    @test nl_function.exp_part[1].convexity == :convex 


    
end