"""
Create the expression graph (DAG)
"""
function create_dag!(model::MOI.AbstractOptimizer)
    model.inner.expression_graph = DAG()
    model.inner.dag_lookup = Dict{Union{Expr, Float64, Int, Symbol}, Tuple{Int, Int}}()
    model.inner.common_sub_expression_dict = Dict{Expr, Vector{Int}}()
    
    if model.inner.is_objective_nl || model.inner.is_objective_quadratic
        id = 0 
        add_to_dag!(model.inner.objective_nl_function, model.inner.expression_graph, 
            model.inner.common_sub_expression_dict, id; dag_lookup = model.inner.dag_lookup)
    end 

    for i in 1:model.inner.num_quadratic_le_constraints 
        nl_function = model.inner.quadratic_nl_function_le[i]
        id = quadratic_le_offset(model) + i
        add_to_dag!(nl_function, model.inner.expression_graph, 
            model.inner.common_sub_expression_dict, id; dag_lookup = model.inner.dag_lookup)
    end 

    for i in 1:model.inner.num_quadratic_ge_constraints 
        nl_function = model.inner.quadratic_nl_function_ge[i]
        id = quadratic_ge_offset(model) + i
        add_to_dag!(nl_function, model.inner.expression_graph, 
            model.inner.common_sub_expression_dict, id; dag_lookup = model.inner.dag_lookup)
    end 

    for i in 1:model.inner.num_quadratic_eq_constraints 
        nl_function = model.inner.quadratic_nl_function_eq[i]
        id = quadratic_eq_offset(model) + i
        add_to_dag!(nl_function, model.inner.expression_graph, 
            model.inner.common_sub_expression_dict, id; dag_lookup = model.inner.dag_lookup)
    end 

    for i in 1:model.inner.num_nlp_constraints
        nl_function = model.inner.nl_function[i]
        id = nlp_offset(model) + i
        add_to_dag!(nl_function, model.inner.expression_graph, 
            model.inner.common_sub_expression_dict, id; dag_lookup = model.inner.dag_lookup)
    end 

    update_max_depth!(model.inner.expression_graph)
    
    return 
end

"""
Update max depth of the DAG 
""" 
function update_max_depth!(dag::DAG)
    depth = [d for d in keys(dag.vertices)]
    (length(depth) == 0) && (return)
    dag.max_depth = maximum(depth)
    return
end 

"""
Add an NLFunction to a DAG 
"""
function add_to_dag!(nl_function::NLFunction, 
    dag::DAG, 
    common_sub_expression_dict::Dict{Expr,Vector{Int}}, 
    constraint_id::Int; 
    dag_lookup = Dict{Union{Expr, Float64, Int, Symbol}, Tuple{Int, Int}}())

    fnames = fieldnames(NLFunction)

    for fname in fnames
        (fname == :linear_part || fname == :constant_part) && (continue)
        f_value = getfield(nl_function, fname)
        if ~isa(f_value, Nothing)
            for i in 1:length(f_value)
                (f_value[i].expression[1] == 0.0) && (continue)
                expr = f_value[i].expression[2]
                depth, id = append_dag!(expr, dag = dag, expr_dict = common_sub_expression_dict, 
                    id = constraint_id, dag_lookup = dag_lookup)
                dag_lookup[expr] = (depth, id)
            end 
        end 
    end 

    return 
end 

"""
Add an nonlinear term, represented by an expression tree to a DAG 
""" 
function append_dag!(expr; 
    dag::DAG = DAG(), 
    expr_dict::Dict{Expr, Vector{Int}} = Dict{Expr,Vector{Int}}(), 
    id::Int = 1, 
    parent::Union{Nothing, Union{Symbol, Expr, Float64, Int}} = nothing, 
    dag_lookup = Dict{Union{Expr, Float64, Int, Symbol}, Tuple{Int, Int}}())::Tuple{Int, Int}

    # handle the case when expr is a number or a single variable
    if isa(expr, Float64) || isa(expr, Int) || 
        (try expr.head == :ref catch y false end) 
        if haskey(dag_lookup, expr)
            depth, position = dag_lookup[expr]
            if ~isa(parent, Nothing)
                push!(dag.vertices[depth][position].parents, parent)
            end
            return dag_lookup[expr]  
        else 
            current_vertex_type = :constant 
            (try expr.head == :ref catch y false end) && (current_vertex_type = :variable)
            current_depth = 0 
            current_vertex = expr
            current_children = Vector{Union{Symbol, Expr, Float64, Int}}() 
            (isa(parent, Nothing)) && (error(LOGGER, "DAG vertex $current_vertex is isolated"))
            current_parents = Vector{Union{Symbol, Expr, Float64, Int}}() 
            push!(current_parents, parent)
            current_dag_vertex = DAGVertex(current_vertex_type, current_depth, current_vertex, current_children, current_parents)
            (~haskey(dag.vertices, 0)) && (dag.vertices[0] = DAGVertex[])
            push!(dag.vertices[current_depth], current_dag_vertex)
            dag_lookup[expr] = (current_depth, length(dag.vertices[current_depth]))
            return current_depth, length(dag.vertices[current_depth])
        end 
    end 

    # common sub-expression check (very important to control the size of the DAG)  
    if haskey(expr_dict, expr)
        if ~haskey(dag_lookup, expr)
            error(LOGGER, "dag_lookup dictionary should contain the expression $expr, since it is contained in expr_dict")
        else 
            depth, position = dag_lookup[expr]
            if ~isa(parent, Nothing)
                push!(dag.vertices[depth][position].parents, parent)
            end
            push!(expr_dict[expr], id)
            return dag_lookup[expr]
        end 
    else 
        expr_dict[expr] = Vector{Int}()
    end 
    push!(expr_dict[expr], id)

    # children depth to compute current DAGVertex depth
    children_depth = Vector{Int}()

    # prepare a DAGVertex for the expression
    current_vertex_type = :operator 
    current_depth = typemax(Int)
    current_vertex = expr
    current_children = Vector{Union{Symbol, Expr, Float64, Int}}() 
    current_parents = Vector{Union{Symbol, Expr, Float64, Int}}() 
    (~isa(parent, Nothing)) && (push!(current_parents, parent))
    current_dag_vertex = DAGVertex(current_vertex_type, current_depth, current_vertex, current_children, current_parents)
    
    # recursive call to the arguments
    for i in 2:length(expr.args)
        push!(current_children, expr.args[i])
        depth, position = append_dag!(expr.args[i], dag = dag, expr_dict = expr_dict, 
            id = id, parent = expr, dag_lookup = dag_lookup)
        push!(children_depth, depth)
    end 

    # update depth and add the vertex to DAG
    current_depth = maximum(children_depth) + 1
    current_dag_vertex.depth = current_depth 
    (~haskey(dag.vertices, current_depth)) && (dag.vertices[current_depth] = DAGVertex[])
    push!(dag.vertices[current_depth], current_dag_vertex)
    dag_lookup[expr] = (current_depth, length(dag.vertices[current_depth]))
    
    return dag_lookup[expr]

end 