"""
Utilities for performing operations on nonlinear expression trees
""" 

"""
Clean the nonlinear expressions
"""
function clean_nl_expressions!(model::MOI.AbstractOptimizer)
    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)
        for expr in model.inner.nl_constraint_expr
            expr_flatten_constant_subtree(expr)
            expr_separate_sign_multilinear(expr)
            expr_aggregate_coeff_multilinear(expr)
        end
        
        if model.inner.is_objective_nl
            expr_flatten_constant_subtree(model.inner.objective_expr)
            expr_separate_sign_multilinear(model.inner.objective_expr)
            expr_aggregate_coeff_multilinear(model.inner.objective_expr)
        end 
    end 

    return
end

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

"""
Break an expression tree into sum of pieces of expression trees;
returns a vector of Tuple{Float64, Union{Expr, Symbol, Float64, Int}}, 
where the first element is the coefficient and the second is Union{Expr, Symbol, Float64, Int}.
This is a core function for all of Alpine.jl
"""
function expr_disaggregate(expr; coeff=1.0)
    expr_tree = Vector{Tuple{Float64, Union{Expr, Symbol, Float64, Int}}}()
    
    if isa(expr, Float64) || isa(expr, Int) 
        push!(expr_tree, (coeff, expr))
        return expr_tree
    end

    if expr.args[1] in [:(<=), :(>=), :(==)]
        append!(expr_tree, expr_disaggregate(expr.args[2], coeff=coeff))
        return expr_tree 
    end
    
    if expr.head == :ref 
        push!(expr_tree, (coeff, expr))
        return expr_tree
    end

    if expr.args[1] == :-
        # (:-) is either be unary or binary operator - both cases handled here
        if length(expr.args) == 2
            append!(expr_tree, expr_disaggregate(expr.args[2], coeff=-1*coeff))
        else 
            append!(expr_tree, expr_disaggregate(expr.args[2], coeff=coeff))
            append!(expr_tree, expr_disaggregate(expr.args[3], coeff=-1*coeff))
        end
    
    elseif expr.args[1] == :+
        # (:+) is n-ary operator - all cases handled here
        for i in 2:length(expr.args)
            append!(expr_tree, expr_disaggregate(expr.args[i], coeff=coeff))
        end

    elseif isa(expr.args[1], Float64) || isa(expr.args[1], Int)
        # numbers handled here
        push!(expr_tree, (coeff, expr.args[1]))
    
    elseif expr.head == :call  
        # :call is handled here; [:+, :-, :*] handled separately 
        # for :*, it is assumed coefficient aggregation is already performed, 
        # otherwise this function will result in faulty disaggregation
        if expr.args[1] in [:+, :-]
            append!(expr_tree, expr_disaggregate(expr, coeff=coeff))
        elseif expr.args[1] == :* && (isa(expr.args[2], Float64) || isa(expr.args[2], Int))
            cleaned_expr = deepcopy(expr)
            constant = cleaned_expr.args[2]
            deleteat!(cleaned_expr.args, 2)
            if length(cleaned_expr.args) == 2
                push!(expr_tree, (coeff * constant, cleaned_expr.args[2]))
            else 
                push!(expr_tree, (coeff * constant, cleaned_expr))
            end 
        else 
            push!(expr_tree, (coeff, expr))
        end

    end 

    return expr_tree

end