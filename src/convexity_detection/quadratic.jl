"""
    get_convexity(quadratic_function::SQF)::Symbol 

Convexity detection for quadratic function (algebraic check)
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
    get_convexity(quadratic_function::SQF, term::SQT)::Symbol 

Checking if a quadratic term is factorable in a quadratic function and deduce convexity (algebraic)
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
    get_square_term(quadratic_function::SQF, id::VI)::Union{SQT, Nothing}

Get square term by variable index in a quadratic function, returns nothing if 
not found
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
    get_convexity(Q::SparseArrays.SparseMatrixCSC{Float64,Int64})::Symbol 

Returns :convex if the minimum eigen value of Q > 0, :concave if maximum eigen 
value of Q < 0, :undet otherwise.
"""
function get_convexity(Q::SparseArrays.SparseMatrixCSC{Float64,Int64})::Symbol 
    eigen_values = eigen(Matrix(Q))
    max_eigen_value = maximum(eigen_values)
    min_eigen_value = minimum(eigen_values)
    if min_eigen_value > 0 
        return :convex 
    elseif max_eigen_value < 0 
        return :concave 
    else 
        return :undet 
    end
end 

"""
    run_quadratic_convexity_detection(model::MOI.AbstractOptimizer)

This function runs the algebriac and PSD test on all the quadratic constraints.
"""
function run_quadratic_convexity_detection!(model::MOI.AbstractOptimizer)
    for i in 1:model.inner.num_quadratic_constraints
        quadratic_function = model.quadratic_constraints[i][1]
        constraint_set = model.quadratic_constraints[i][2]
        Q = model.inner.quadratic_matrix[i].Q 
        is_function_convex = get_convexity(quadratic_function)
        if is_function_convex == :undet 
            is_function_convex = get_convexity(Q)
        end
        model.inner.quadratic_function_convexity[i] = is_function_convex
        if is_function_convex == :convex
            if isa(constraint_set, MOI.LessThan{Float64}) 
                model.inner.quadratic_constraint_convexity[i] = :convex 
            elseif isa(constraint_set, MOI.GreaterThan{Float64})
                model.inner.quadratic_constraint_convexity[i] = :concave 
            end
        elseif is_function_convex == :concave 
            if isa(constraint_set, MOI.LessThan{Float64}) 
                model.inner.quadratic_constraint_convexity[i] = :concave
            elseif isa(constraint_set, MOI.GreaterThan{Float64})
                model.inner.quadratic_constraint_convexity[i] = :convex
            end
        end
    end 

    if model.sense != MOI.FEASIBILITY_SENSE && model.inner.is_objective_quadratic
        is_function_convex = get_convexity(model.objective)
        Q = model.inner.quadratic_matrix_objective.Q 
        if is_function_convex == :undet 
            is_function_convex = get_convexity(Q)
        end
        model.inner.objective_function_convexity = is_function_convex
        if is_function_convex == :convex && model.sense == MOI.MIN_SENSE 
            model.inner.objective_convexity = :convex
        end 
        if is_function_convex == :concave && model.sense == MOI.MAX_SENSE 
            model.inner.objective_convexity = :convex
        end 
    end

    return
end 

"""
    run_quadratic_nl_convexity_detection(model::MOI.AbstractOptimizer)

This function runs the PSD test on the quadratic parts of the nonlinear constraints. 
"""
function run_quadratic_nl_convexity_detection!(model::MOI.AbstractOptimizer)
    for i in 1:model.inner.num_nl_constraints 
        nl_function = model.inner.constraint_nl_function[i]
        quadratic_part = nl_function.quadratic_part 
        Q = model.inner.quadratic_matrix_nl[i]
        if isa(Q, Nothing) || isa(quadratic_part, Nothing)
            continue 
        end
        model.inner.nl_quadratic_matrix_convexity = get_convexity(Q)
    end

    if model.inner.is_objective_nl 
        nl_function = model.inner.objective_nl_function 
        quadratic_part = nl_function.quadratic_part 
        Q = model.inner.quadratic_matrix_objective
        if isa(Q, Nothing) || isa(quadratic_part, Nothing)
            return 
        end
        model.inner.objective_quadratic_matrix_convexity = get_convexity(Q)
    end 
    
    return 
end 