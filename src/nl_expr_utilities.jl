"""
Utilities for performing operations on nonlinear expression trees
""" 

"""
Check if a sub-tree is totally composed of constants
"""
function expr_is_constant(expr)
	(isa(expr, Float64) || isa(expr, Int) || isa(expr, Symbol)) && return true
    is_constant = true

    for i in 1:length(expr.args)
		if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
			continue
		elseif expr.args[i].head == :call
			is_constant = is_constant * expr_is_constant(expr.args[i])
		elseif expr.args[i].head == :ref
			return false
		end
	end

	return is_constant
end

"""
Flatten sub-trees that contain only constant values
"""
function expr_flatten_constant_subtree(expr)    
    for i in 1:length(expr.args)
		if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
			continue
		elseif expr.args[i].head == :call
			(expr_is_constant(expr.args[i])) && (expr.args[i] = eval(expr.args[i]))
			if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
				continue
			else
				expr_flatten_constant_subtree(expr.args[i])
			end
		end
	end
    
    return
end 

"""
Multilinear term negative sign separation 
"""
function expr_separate_sign_multilinear(expr)
    for i in 2:length(expr.args)
        if !isa(expr.args[i], Float64) && !isa(expr.args[i], Int) 
            if length(expr.args[i].args) == 2 && expr.args[i].head == :call 
                if expr.args[1] == :*
                    if expr.args[i].args[1] == :- 
                        push!(expr.args, -1.0)
                        expr.args[i] = expr.args[i].args[2]
                    elseif expr.args[i].args[1] == :+
                        expr.args[i] = expr.args[i].args[2]
                    end 
                end 
            elseif expr.args[i].head == :call 
                expr_separate_sign_multilinear(expr.args[i])
            end 
        end 
    end 

    return
end

"""
Multilinear term coefficient aggregration 
"""
function expr_aggregate_coeff_multilinear(expr)
    if expr.args[1] == :* 
        constant = 1.0
        indexes = []
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                constant *= expr.args[i]
                push!(indexes, i)
            end 
        end 
        deleteat!(expr.args, indexes)
        (constant != 1.0) && insert!(expr.args, 2, constant) 
    end 

    for i in 1:length(expr.args)
        if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
			continue
        elseif expr.args[i].head == :call
            expr_aggregate_coeff_multilinear(expr.args[i])
        end 
    end 

    return 
end 