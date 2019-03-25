"""
This file contains all the functions to classify the type of an expression tree 
"""

"""
Check if term is quadratic 
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
Check if term is x^a, where a != 2
"""
function expr_is_power(expr)::Bool 
    (length(expr.args) > 3 || length(expr.args) < 2) && (return false)
    (expr_is_quadratic(expr)) && (return false)
    is_power = (expr.args[1] == :^ || expr.args[1] == :sqrt)

    if is_power && expr.args[1] == :^ 
        is_power = is_power && 
        try expr.args[2].head == :ref && 
            (isa(expr.args[3], Float64) || isa(expr.args[3], Int))
        catch y 
            false 
        end 
    end 

    if is_power && expr.args[1] == :sqrt 
        is_power = is_power && 
        try expr.args[2].head == :ref && length(expr.args) == 2
        catch y 
            false 
        end 
    end 

    return is_power 
end

"""
Check if term is bilinear 
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
Check if term is multilinear
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

"""
Check if term is polynomial
"""
function expr_is_polynomial(expr)::Bool 
    (expr.args[1] != :*) && (return false)
    is_polynomial = true
    variable_powers = Dict{Int, Int}()

    for i in 2:length(expr.args)
        
    end 

    return is_polynomial

end 

"""
Check if term is signomial 
""" 
function expr_is_signomial(expr)::Bool 
    (expr.args[1] != :*) && (return false)
    is_signomial = true

    return is_signomial
end 

"""
Check if term is abs(x)
"""
function expr_is_abs(expr)::Bool 
    (expr.args[1] != :abs) && (return false)
    is_abs = true

    return is_abs 
end 

"""
Check if term is sin(x), cos(x), tan(x)
"""
function expr_is_trignometric(expr)::Bool 
    operations = [:sin, :cos, :tan]
    (~(expr.args[1] in operations)) && (return false)
    is_trigonometric = true

    return is_trigonometric 
end 

"""
Check if term is log(x)
"""
function expr_is_log(expr)::Bool 
    (expr.args[1] != :log) && (return false)
    is_log = true
    
    return is_log 
end 

"""
Check if term is exp(x)
"""
function expr_is_exp(expr)::Bool 
    (expr.args[1] != :exp) && (return false)
    is_exp = true
    
    return is_exp
end 