
@testset "DAG construction tests" begin
    
    printstyled("\nTesting DAG construction ...\n", color=:blue, bold=true)
    
    solver_options = DefaultTestSolver()
    m = Model(with_optimizer(
        Alpine.Optimizer, solver_options)
    )

    @variable(m, x[1:3])

    dag = Alpine.DAG()
    dag_lookup = Dict{Union{Int, Symbol, Float64, Expr}, Tuple{Int, Int}}()
    expr_dict = Dict{Expr,Vector{Int}}()

    expressions = [
        :(log(x[1]^2 + x[2]^2) * x[3]), # 1
        :(x[1]^2 * x[2]^2), # 2
        :(sin(x[1]) * x[2]), # 3
        :(x[3]^1.5 * abs(x[3])), # 4 
        :(x[3]^1.5 * x[1] * x[2]), # 5
        :(x[1] * x[2]), # 6
        :(abs(log(x[1] * x[2]))) # 7
    ]

    for i in 1:length(expressions)
        expr = expressions[i]
        Alpine.append_dag!(expr, dag = dag, expr_dict = expr_dict, id = i, parent = nothing, dag_lookup = dag_lookup)
    end
    Alpine.update_max_depth!(dag)
    
    @test length(expr_dict) == 15
    @test expr_dict[:(x[3] ^ 1.5)] == [4, 5]
    @test expr_dict[:(x[1] * x[2])] == [6, 7]
    @test dag.max_depth == 4
    @test length(dag.vertices[0]) == 5
    @test length(dag.vertices[1]) == 6
    @test length(dag.vertices[2]) == 6
    @test length(dag.vertices[3]) == 2
    @test length(dag.vertices[4]) == 1
    @test dag_lookup[:(x[1] ^ 2 * x[2] ^ 2)][1] == 2
    @test dag_lookup[:(x[1] ^ 2 * x[2] ^ 2)][2] == 2
    @test length(dag.vertices[2][2].parents) == 0

end 