"""
    populate_nl_convexity!(model::MOI.AbstractOptimizer)

This function determines the convexity of the nonlinear function and the
nonlinear constraint as a whole and updates these fields
"""

function populate_nl_convexity!(model::MOI.AbstractOptimizer)
    dag_lookup = model.inner.dag_lookup
    dag = model.inner.expression_graph
    for i in 1:model.inner.num_nl_constraints 
        nl_function = model.inner.nl_function[i]
        disaggregated_convexity = Symbol[]
        if ~isa(nl_function.quadratic_part, Nothing)
            push!(disaggregated_convexity, model.inner.nl_quadratic_matrix_convexity[i])
        end
        lb = model.nlp_data.constraint_bounds[i].lower 
        ub = model.nlp_data.constraint_bounds[i].upper 
        constraint_type = :NaN 
        if lb == ub
            constraint_type = :(==)
        elseif !isinf(lb) 
            constraint_type = :(>=)
        else 
            constraint_type = :(<=)
        end
        fnames = fieldnames(NLFunction)
        for fname in fnames 
            f_value = getfield(nl_function, fname)
            (isa(f_value, Nothing)) && (continue)
            (fname in [:constant_part, :linear_part, :quadratic_part]) && (continue)
            for alpine_expression in f_value 
                coeff, expr = alpine_expression.expression
                depth, position = dag_lookup[expr]
                dag_vertex = dag.vertices[depth][position] 
                if convex(dag_vertex.convexity) && coeff >= 0.0
                    alpine_expression.convexity = :convex 
                    push!(disaggregated_convexity, :convex)
                elseif convex(dag_vertex.convexity) && coeff <= 0.0 
                    alpine_expression.convexity = :concave 
                    push!(disaggregated_convexity, :concave)
                elseif concave(dag_vertex.convexity) && coeff >= 0.0 
                    alpine_expression.convexity = :concave 
                    push!(disaggregated_convexity, :concave)
                elseif concave(dag_vertex.convexity) && coeff <= 0.0 
                    alpine_expression.convexity = :convex
                    push!(disaggregated_convexity, :convex)
                else 
                    push!(disaggregated_convexity, :undet)
                end
            end 
        end 
        is_convex = convex(disaggregated_convexity)
        is_concave = concave(disaggregated_convexity)

        nl_function_convexity = get_convexity(is_convex, is_concave)
        model.inner.nl_function_convexity[i] = nl_function_convexity

        if convex(nl_function_convexity) && constraint_type == :(<=)
            model.inner.nl_constraint_convexity[i] = :convex 
        elseif concave(nl_function_convexity) && constraint_type == :(>=)
            model.inner.nl_constraint_convexity[i] = :convex 
        elseif convex(nl_function_convexity) && constraint_type == :(>=) 
            model.inner.nl_constraint_convexity[i] = :concave 
        elseif concave(nl_function_convexity) && constraint_type == :(<=)
            model.inner.nl_constraint_convexity[i] = :concave 
        else 
            model.inner.nl_constraint_convexity[i] = :undet
        end
    end 

    return
end 

"""
    populate_objective_convexity!(model::MOI.AbstractOptimizer)

This function determines the convexity of the nonlinear objective function 
"""
function populate_objective_convexity!(model::MOI.AbstractOptimizer)
    dag_lookup = model.inner.dag_lookup
    dag = model.inner.expression_graph
    if model.sense != MOI.FEASIBILITY_SENSE && model.inner.is_objective_nl 
        nl_function = model.inner.objective_nl_function
        disaggregated_convexity = Symbol[]
        if ~isa(nl_function.quadratic_part, Nothing)
            push!(disaggregated_convexity, model.inner.objective_quadratic_matrix_convexity)
        end
        
        fnames = fieldnames(NLFunction)
        for fname in fnames 
            f_value = getfield(nl_function, fname)
            (isa(f_value, Nothing)) && (continue)
            (fname in [:constant_part, :linear_part, :quadratic_part]) && (continue)
            for alpine_expression in f_value 
                coeff, expr = alpine_expression.expression
                depth, position = dag_lookup[expr]
                dag_vertex = dag.vertices[depth][position] 
                if convex(dag_vertex.convexity) && coeff >= 0.0
                    alpine_expression.convexity = :convex 
                    push!(disaggregated_convexity, :convex)
                elseif convex(dag_vertex.convexity) && coeff <= 0.0 
                    alpine_expression.convexity = :concave 
                    push!(disaggregated_convexity, :concave)
                elseif concave(dag_vertex.convexity) && coeff >= 0.0 
                    alpine_expression.convexity = :concave 
                    push!(disaggregated_convexity, :concave)
                elseif concave(dag_vertex.convexity) && coeff <= 0.0 
                    alpine_expression.convexity = :convex
                    push!(disaggregated_convexity, :convex)
                else 
                    push!(disaggregated_convexity, :undet)
                end
            end 
        end 
        is_convex = convex(disaggregated_convexity)
        is_concave = concave(disaggregated_convexity)

        nl_function_convexity = get_convexity(is_convex, is_concave)
        model.inner.objective_function_convexity = nl_function_convexity

        if convex(nl_function_convexity) && model.sense == MOI.MIN_SENSE 
            model.inner.objective_convexity = :convex
        elseif concave(nl_function_convexity) && model.sense == MOI.MAX_SENSE 
            model.inner.objective_convexity = :convex
        else
            model.inner.objective_convexity = :undet
        end 
    end 

    return

end 