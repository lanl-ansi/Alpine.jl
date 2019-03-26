"""
Convexity detection for quadratic functions 
"""
function get_convexity(quadratic_function::SQF)::Symbol 

    num_terms = length(quadratic_function.quadratic_terms)
    full_convexity = get_convexity(quadratic_function, quadratic_function.quadratic_terms[1])
    (full_convexity == :undet) && (return :undet)

    for i in 2:num_terms 
        term_convexity = get_convexity(quadratic_function, quadratic_function.quadratic_terms[i])
        (term_convexity == :undet) && (return :undet)
        (full_convexity == :convex && term_convexity == :concave) && (return :undet)
        (full_convexity == :concave && term_convexity == :convex) && (return :undet)
    end 

    return full_convexity
end 

"""
Checking if a quadratic term is factorable in a quadratic function and deduce convexity
"""
function get_convexity(quadratic_function::SQF, term::SQT)::Symbol
    
    if term.variable_index_1 == term.variable_index_2 
        if term.coefficient >= 0.0
            return :convex 
        else 
            return :concave
        end 
    end 

    term_sqr_1 = get_square_term(quadratic_function, term.variable_index_1)
    term_sqr_2 = get_square_term(quadratic_function, term.variable_index_2)

    if term_sqr_1 != nothing && term_sqr_2 != nothing 
        coeff_1 = term_sqr_1.coefficient
        coeff_2 = term_sqr_2.coefficient
        if sign(coeff_1) == sign(coeff_2)
            condition = (2 * sqrt(coeff_1 * coeff_2) >= term.coefficient)
            if condition 
                if sign(coeff_1) == 1
                    return :convex 
                else 
                    return :concave 
                end 
            end 
        end
    end 

    return :undet

end 

"""
Get square term by variable index in a quadratic function 
"""
function get_square_term(quadratic_function::SQF, id::VI)::Union{SQT, Nothing}
    for term in quadratic_function.quadratic_terms
        if term.variable_index_1 == id && term.variable_index_2 == id 
            return term 
        end
    end 
    return nothing
end 


"""
Convexity detection for quadratic objective and constraints  
"""

function run_convexity_detection!(model::MOI.AbstractOptimizer)
    for i in 1:model.inner.num_quadratic_le_constraints
        quadratic_function = model.quadratic_le_constraints[i][1]
        constraint_set = model.quadratic_le_constraints[i][2]
        is_function_convex = get_convexity(quadratic_function)
        model.inner.quadratic_function_convexity[:quadratic_le][i] = is_function_convex
        if is_function_convex == :convex
            model.inner.quadratic_constraint_convexity[:quadratic_le][i] = :convex 
        elseif is_function_convex == :concave 
            model.inner.quadratic_constraint_convexity[:quadratic_le][i] = :concave 
        else 
            model.inner.quadratic_constraint_convexity[:quadratic_le][i] = :undet
        end
    end 

    for i in 1:model.inner.num_quadratic_ge_constraints
        quadratic_function = model.quadratic_ge_constraints[i][1]
        constraint_set = model.quadratic_ge_constraints[i][2]
        is_function_convex = get_convexity(quadratic_function)
        model.inner.quadratic_function_convexity[:quadratic_ge][i] = is_function_convex
        if is_function_convex == :convex
            model.inner.quadratic_constraint_convexity[:quadratic_ge][i] = :concave
        elseif is_function_convex == :concave 
            model.inner.quadratic_constraint_convexity[:quadratic_ge][i] = :convex 
        else 
            model.inner.quadratic_constraint_convexity[:quadratic_ge][i] = :undet
        end
    end 

    for i in 1:model.inner.num_quadratic_eq_constraints
        quadratic_function = model.quadratic_eq_constraints[i][1]
        is_function_convex = get_convexity(quadratic_function)
        model.inner.quadratic_function_convexity[:quadratic_eq][i] = is_function_convex
    end 

    if model.sense != MOI.FEASIBILITY_SENSE && model.inner.is_objective_quadratic
        is_function_convex = get_convexity(model.objective)
        if is_function_convex == :convex && model.sense == MOI.MIN_SENSE 
            model.inner.objective_convexity = :convex
        end 
        if is_function_convex == :concave && model.sense == MOI.MAX_SENSE 
            model.inner.objective_convexity = :convex
        end 
    end

end 
