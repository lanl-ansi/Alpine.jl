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
function expr_disaggregate(expr; coeff=1.0)::Vector{Tuple{Int, Union{Expr, Symbol, Float64, Int}}}
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

"""
Create NLFunction() from disaggregated expression 
""" 
function create_nl_function(disaggregated_expr::Vector{Tuple{Int, Union{Expr, Symbol, Float64, Int}}})::NLFunction 

    nl_function = NLFunction()
    linear_part = AlpineExpr[]
    quadratic_part = AlpineExpr[]    
    power_part = AlpineExpr[]            
    bilinear_part = AlpineExpr[]
    multilinear_part = AlpineExpr[]
    polynomial_part = AlpineExpr[]
    signomial_part = AlpineExpr[]            
    abs_part = AlpineExpr[]
    trigonometric_part = AlpineExpr[]
    log_part = AlpineExpr[]
    exp_part = AlpineExpr[]
    other_part = AlpineExpr[]
    constant_part = AlpineExpr[]

    linear_expr_dict = Dict()
    quadratic_expr_dict = Dict()
    power_expr_dict = Dict()
    bilinear_expr_dict = Dict() 
    multilinear_expr_dict = Dict() 
    polynomial_expr_dict = Dict() 
    signomial_expr_dict = Dict() 
    abs_expr_dict = Dict()
    trigonometric_expr_dict = Dict()
    log_expr_dict = Dict() 
    exp_expr_dict = Dict() 
    other_part_dict = Dict()

    for (coeff, expr) in disaggregated_expr 
        
        # constant term
        if isa(expr, Float64) || isa(expr, Int)
            push!(constant_part, AlpineExpr((coeff, expr), :convex))
            
        # linear term
        elseif expr.head == :ref 
            if haskey(linear_expr_dict, expr)
                index = linear_expr_dict[expr]
                new_coeff = linear_part[index].expression[1] + coeff 
                linear_part[index] = AlpineExpr((new_coeff, expr), :convex)
            else 
                push!(linear_part, AlpineExpr((coeff, expr), :convex))
                linear_expr_dict[expr] = length(linear_part)
            end 
        
        # quadratic term
        elseif expr_is_quadratic(expr)
            if haskey(quadratic_expr_dict, expr)
                index = quadratic_expr_dict[expr]
                new_coeff = quadratic_part[index].expression[1] + coeff
                convexity = (new_coeff >= 0) ? :convex : :concave 
                quadratic_part[index] = AlpineExpr((new_coeff, expr), convexity)
            else 
                convexity = (coeff >= 0) ? :convex : :concave
                push!(quadratic_part, AlpineExpr((coeff, expr), convexity))
                quadratic_expr_dict[expr] = length(quadratic_part)
            end 
        
        # bilinear term 
        elseif expr_is_bilinear(expr)
            if haskey(bilinear_expr_dict, expr)
                index = bilinear_expr_dict[expr]
                new_coeff = bilinear_part[index].expression[1] + coeff 
                bilinear_part[index] = AlpineExpr((new_coeff, expr), :undet)
            else 
                push!(bilinear_part, AlpineExpr((coeff, expr), :undet))
                bilinear_expr_dict[expr] = length(bilinear_part)
            end 
        
        # multilinear term
        elseif expr_is_multilinear(expr)
            if haskey(multilinear_expr_dict, expr)
                index = multilinear_expr_dict[expr]
                new_coeff = multilinear_part[index].expression[1] + coeff 
                multilinear_part[index] = AlpineExpr((new_coeff, expr), :undet)
            else 
                push!(multilinear_part, AlpineExpr((coeff, expr), :undet))
                multilinear_expr_dict[expr] = length(multilinear_part)
            end 
        
        # other terms
        else 
            if haskey(other_expr_dict, expr) 
                index = other_expr_dict[expr]
                new_coeff = other_part[index].expression[1] + coeff 
                other_part[index] = AlpineExpr((new_coeff, expr), :undet)
            else 
                push!(other_part, AlpineExpr((coeff, expr), :undet))
                other_expr_dict[expr] = length(other_part)
            end 

        end  
    end 

    nl_function.constant_part = constant_part
    nl_function.linear_part = linear_part 
    nl_function.quadratic_part = quadratic_part
    nl_function.bilinear_part = bilinear_part
    nl_function.multilinear_part = multilinear_part
    
    return nl_function
end
