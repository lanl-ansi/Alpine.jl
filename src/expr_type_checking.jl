"""
This file contains all the functions to classify the type of an expression tree 
"""

"""
Check if expression is quadratic 
"""
function expr_is_quadratic(expr)::Bool 
    (length(expr.args) != 3) && (return false)
    is_quadratic = (expr.args[1] == :^ && expr.args[3] == 2.0) 
    is_quadratic = is_quadratic || (expr.args[1] == :* && expr.args[2] == expr.args[3])
    
    if is_quadratic && expr.args[1] == :* 
        new_expr = :(expr.args[2]^2.0)
        expr = new_expr 
    end 

    return is_quadratic 
end 

"""
Check if expression is bilinear 
"""
function expr_is_bilinear(expr)::Bool 
    (length(expr.args) != 3) && (return false)
    (expr_is_quadratic(expr)) && (return false)
    (expr.args[1] != :*) && (return false)
    is_bilinear = true
    for i in 2:3
        is_bilinear = is_bilinear && 
        try expr.args[i].head == :ref 
        catch y 
            false 
        end 
    end
    
    if is_bilinear 
        index_1 = expr.args[2].args[2].value
        index_2 = expr.args[3].args[2].value
        expr.args[2].args[2] = VI(min(index_1, index_2))
        expr.args[3].args[2] = VI(max(index_1, index_2))
    end 
    
    return is_bilinear 
end

"""
Check if expression is multilinear
"""
function expr_is_multilinear(expr)::Bool 
    (expr.args[1] != :*) && (return false)
    is_multilinear = true
    variable_indexes = Int[]
    for i in 2:length(expr.args)
        is_multilinear = is_multilinear && 
        try expr.args[i].head == :ref 
        catch y 
            false 
        end
        (is_multilinear) && (push!(variable_indexes, expr.args[i].args[2].value))
        (~is_multilinear) && (return false)
    end 

    sort!(variable_indexes)

    if length(variable_indexes) != length(unique(variable_indexes)) 
        return false
    end 

    if is_multilinear 
        for i in 2:length(expr.args)
            expr.args[i].args[2] = VI(variable_indexes[i-1])
        end 
    end 

    return is_multilinear

end 