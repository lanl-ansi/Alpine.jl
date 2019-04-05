"""
Generic utilities used in Alpine.jl 
"""
function matrix_from_quadratic_terms(terms::Vector{SQT})
    variable_to_index_map = Dict{VI, Int}()
    index_to_variable_map = Dict{Int, VI}()
    n = 0
    for term in terms
        for variable in (term.variable_index_1, term.variable_index_2)
            if !(variable in keys(variable_to_index_map))
                n += 1
                variable_to_index_map[variable] = n
                index_to_variable_map[n] = variable
            end
        end
    end

    I = Int[]
    J = Int[]
    V = Float64[]
    
    for term in terms
        i = variable_to_index_map[term.variable_index_1]
        j = variable_to_index_map[term.variable_index_2]
        push!(I, i)
        push!(J, j)
        push!(V, term.coefficient)
        if i != j
            push!(I, j)
            push!(J, i)
            push!(V, term.coefficient)
        end
    end

    # Duplicate terms are summed together in `sparse`
    Q = sparse(I, J, V, n, n)
    
    return Q, index_to_variable_map
end

function matrix_from_nl_function(nl_function::NLFunction)
    quadratic_part = nl_function.quadratic_part 
    bilinear_part = nl_function.bilinear_part 
    if isa(quadratic_part, Nothing) && isa(bilinear_part, Nothing)
        return nothing, nothing 
    end
    variable_to_index_map = Dict{VI, Int}()
    index_to_variable_map = Dict{Int, VI}()
    n = 0
    if ~isa(quadratic_part, Nothing)
        for quad_term in quadratic_part 
            coeff, expr = quad_term.expression 
            variable = expr.args[2].args[2]
            if !(variable in keys(variable_to_index_map))
                n += 1
                variable_to_index_map[variable] = n
                index_to_variable_map[n] = variable
            end
        end 
    end 

    if ~isa(bilinear_part, Nothing)
        for bilinear_term in bilinear_part 
            coeff, expr = bilinear_term.expression 
            variable_1 = expr.args[2].args[2]
            variable_2 = expr.args[3].args[2]
            for variable in (variable_1, variable_2)
                if !(variable in keys(variable_to_index_map))
                    n += 1
                    variable_to_index_map[variable] = n
                    index_to_variable_map[n] = variable
                end
            end
        end
    end 

    I = Int[]
    J = Int[]
    V = Float64[]

    if ~isa(quadratic_part, Nothing)
        for quad_term in quadratic_part 
            coeff, expr = quad_term.expression 
            variable = expr.args[2].args[2]
            i = variable_to_index_map[variable]
            push!(I, i)
            push!(J, i)
            push!(V, coeff * 2.0)
        end 
    end 

    if ~isa(bilinear_part, Nothing)
        for bilinear_term in bilinear_part 
            coeff, expr = bilinear_term.expression 
            variable_1 = expr.args[2].args[2]
            variable_2 = expr.args[3].args[2]
            i = variable_to_index_map[variable_1]
            j = variable_to_index_map[variable_2]
            push!(I, i)
            push!(J, j)
            push!(V, coeff)
            push!(I, j)
            push!(J, i)
            push!(V, coeff)
        end
    end 

    # Duplicate terms are summed together in `sparse`
    Q = sparse(I, J, V, n, n)
    
    return Q, index_to_variable_map
end 
