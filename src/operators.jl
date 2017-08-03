#  Proposed built-in structure patterns
#   Plan for the following non-linear subtree structure
#   1. bilinear tree
#           (:*)        x*y
#           /  \
#         y     x
#   2. monomial tree
#           (:^)        x^2
#           /  \
#         x     2
#   3. power tree       x^a, a>2
#           (:^)
#           /  \
#         x    >2
#   3. multilinear tree  x*z*y
#           (:*)
#          /| ..\
#        x  z   y
#   4. hierarchical bilinear tree (upon 1)  (x*y)*z
#           (:*)
#           /  \
#         (:*)  z
#         /  \
#        x    y
#   5. sin tree  (can fit 1, 2, 3)   sin(x)
#          (:sin)
#             |
#             x
#   6. cos tree (can fit 1, 2, 3)    cos(x)
#           (:cos)
#             |
#             x
#   7. user-defined tree (can over-ride built-in structural trees)
#
#

"""
    process_expr(expr; kwargs...)

High-level warpper for processing expression with sub-tree operators
"""
function expr_batch_process(m::PODNonlinearModel;kwargs...)

	# 0 : deepcopy data into mip lifted expr place holders
	m.bounding_obj_expr_mip = deepcopy(m.obj_expr_orig)
	for i in 1:m.num_constr_orig
        push!(m.bounding_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
    end

    # 1 : pre-process the negative sign in expressions
    expr_resolve_const(m.bounding_obj_expr_mip)
    expr_resolve_sign(m.bounding_obj_expr_mip)
    expr_flatten(m.bounding_obj_expr_mip)
    for i in 1:m.num_constr_orig
        expr_resolve_const(m.bounding_constr_expr_mip[i])
        expr_resolve_sign(m.bounding_constr_expr_mip[i])
        expr_flatten(m.bounding_constr_expr_mip[i].args[2])
    end

    # 2 : most important
    m.bounding_obj_expr_mip = expr_term_parsing(m.bounding_obj_expr_mip, m)
    for i in 1:m.num_constr_orig
        # Expression parsing and recognizing
        m.bounding_constr_expr_mip[i] = expr_term_parsing(m.bounding_constr_expr_mip[i], m)
    end

    # 3 : extract some information
    for i in keys(m.nonlinear_terms)
        for var in i
            @assert isa(var.args[2], Int)
            if !(var.args[2] in m.all_nonlinear_vars)
                push!(m.all_nonlinear_vars, var.args[2])
            end
        end
    end
    m.all_nonlinear_vars = sort(m.all_nonlinear_vars)
    m.num_var_lifted_mip = length(m.nonlinear_terms)

    # 4 : convert expression graph into affine functions
    populate_bounding_mip_oc(m)

    return m
end

"""
    expr_constr_parsing(expr, m::PODNonlinearModel)

Recognize structural constraints.
"""
function expr_constr_parsing(expr, m::PODNonlinearModel, idx::Int; kwargs...)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.constr_patterns)
        skip, expr = eval(m.constr_patterns[i])(expr, m)
        skip && return expr
    end

    # Recognize linear constriant
    linear && (m.structural_constr[idx] = :linear)

    # Recognize built-in special structural patterns
    convex, expr = expr_isconvex_constr(expr)
    convex && (m.structural_constr[idx] = :convex)

    # The rest of the constraints will be recognized as :nonlinear
    m.structural_constr[idx] = :nonlinear

    return expr
end

"""
	Check if a sub-tree(:call) is totally composed of constant values
"""
function expr_isconst(expr)

	(isa(expr, Float64) || isa(expr, Int) || isa(expr, Symbol)) && return true

	const_tree = true
	for i in 1:length(expr.args)
		if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
			continue
		elseif expr.args[i].head == :call
			const_tree *= expr_isconst(expr.args[i])
		elseif expr.args[i].head == :ref
			return false
		end
	end

	return const_tree
end

"""
    expr_term_parsing(expr, m::PODNonlinearModel, level=0)

Recognize and process nonlinear terms in an expression
"""
function expr_term_parsing(expr, m::PODNonlinearModel, level=0;options...)

    cnt = 0
    for node in expr.args
        cnt += 1
        if isa(node, Float64) || isa(node, Int) || isa(node, Symbol)
            continue
        elseif node.head == :call
            expr.args[cnt] = expr_term_parsing(node, m, level+1)
        elseif node.head == :ref
            continue
        else
            error("Type issue during expression parsing. ")
        end
    end

    return expr_resolve_term_pattern(expr, m)
end

"""
    expr_resolve_term_pattern(expr, m::PODNonlinearModel)

This function recognizes, stores, and replaces a sub-tree `expr` with available
user-defined/built-in structures patterns. The procedure is creates the required number
of lifted variables based on the patterns that it it trying to recognize.
Then, go through all built-in structures and perform operatins to convexify the problem.

Available structures patterns are:
    * bilinear
    * monomial
    * multi-linears
    * sin
    * cos
    * user-defined

Specific structure pattern information will be described formally.
"""
function expr_resolve_term_pattern(expr, m::PODNonlinearModel; kwargs...)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.term_patterns)
        skip, expr = eval(m.term_patterns[i])(expr, m)
        skip && return expr
    end

    # Recognize all built-in structural patterns
    skip, expr = resolve_bilinear(expr, m)
    skip && return expr

    skip, expr = resolve_monomial(expr, m)
    skip && return expr

    skip, expr = resolve_multilinear(expr, m)
    skip && return expr

    # resolve_sin(expr,m) && return expr
    # resolve_cos(expr,m) && return expr

    return expr # if no structure is detected, simply return the original tree
end

"""
    TODO: docstring
"""
function resolve_bilinear(expr, m::PODNonlinearModel)

    @assert expr.head == :call

    function store_bilinear()
        y_idx = m.num_var_orig + length(keys(m.nonlinear_terms)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[2])))
        m.nonlinear_terms[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                        :id => length(keys(m.nonlinear_terms)) + 1,
                                        :ref => term_key,
                                        :lifted_constr_ref => lifted_constr_ref,
                                        :nonlinear_type => :bilinear,
                                        :convexified => false)
    end

    function lift_bilinear()
        if scalar == 1
            return m.nonlinear_terms[term_key][:lifted_var_ref]
        else
            return Expr(:call, :*, m.nonlinear_terms[term_key][:lifted_var_ref], scalar)
        end
    end

    # Main body
    if (expr.args[1] == :*)  # confirm head (:*)
        # ----- Pattern : coefficient * x * y  ------ #
        # Collect children information for checking
        scalar = 1.0
        var_idxs = []
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false, expr
        end
        # Cofirm detection of patter A and perform store & lifting procedures
        if length(var_idxs) == 2
            (m.log_level) > 99 && println("found bilinear term $expr")
            term_key = [Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[2])]
            if (term_key in keys(m.nonlinear_terms) || reverse(term_key) in keys(m.nonlinear_terms))
                (term_key in keys(m.nonlinear_terms)) ? term_key = term_key : term_key = reverse(term_key)
                return true, lift_bilinear()
            else
                store_bilinear()
                return true, lift_bilinear()
            end
        end
    end

    return false, expr
end

"""
    TODO: docstring
"""
function resolve_multilinear(expr, m::PODNonlinearModel)

    function store_multilinear()
        y_idx = m.num_var_orig + length(keys(m.nonlinear_terms)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        constr_block = "x[$(y_idx)]=="
        for j in 1:length(var_idxs)
            constr_block = string(constr_block, "x[$(var_idxs[j])]")
            if j < length(var_idxs)
                constr_block=string(constr_block, "*")
            end
        end
        lifted_constr_ref = parse(constr_block)
        m.nonlinear_terms[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :id => length(keys(m.nonlinear_terms)) + 1,
                                            :ref => term_key,
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :nonlinear_type => :multilinear,
                                            :convexified => false)
    end

    function lift_multilinear()
        if scalar == 1
            return m.nonlinear_terms[term_key][:lifted_var_ref]
        else
            return Expr(:call, :*, m.nonlinear_terms[term_key][:lifted_var_ref], scalar)
        end
    end

    @assert expr.head == :call
    if (expr.args[1] == :*)     # confirm head (:*)
        # Pattern: coefficients * x * y * z ...
        var_idxs = []
        scalar = 1.0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false, expr
        end
        if length(var_idxs) > 2
            (m.log_level > 99) && println("found multilinear term $expr")
            term_key = Set()
            [push!(term_key, Expr(:ref, :x, idx)) for idx in var_idxs]
            if term_key in keys(m.nonlinear_terms)
                return true, lift_multilinear()
            else
                store_multilinear()
                return true, lift_multilinear()
            end
        end
    end

    return false, expr
end

"""
    TODO: docstring
"""
function resolve_monomial(expr, m::PODNonlinearModel)

    function store_monomial()
        y_idx = m.num_var_orig + length(keys(m.nonlinear_terms)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[1])))
        m.nonlinear_terms[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :id => length(keys(m.nonlinear_terms)) + 1,
                                            :ref => term_key,
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :nonlinear_type => :monomial,
                                            :convexified => false)
    end

    function lift_monomial()
        new_expr = m.nonlinear_terms[term_key][:lifted_var_ref] #TODO: check if this will actually be replaced
        return new_expr
    end

    if (expr.args[1] == :^) && length(expr.args) == 3
        # Pattern: (x)^(2)
        var_idxs = []
        power_scalar = 0
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                power_scalar += expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false, expr
        end
        if length(var_idxs) == 1 && power_scalar == 2.0
            (m.log_level > 99) && println("found monomial term $expr")
            term_key = []
            for i in 1:2
                push!(term_key, Expr(:ref, :x, var_idxs[1]))
            end

            if term_key in keys(m.nonlinear_terms)
                return true, lift_monomial()
            else
                store_monomial()
                return true, lift_monomial()
            end
        end
    end

    return false, expr
end

"""
    TODO: doc
"""
function resolve_sine(expr, m::PODNonlinearModel)

    function store_sin()
    end

    function lift_sin()
    end

    # @assert expr.head == :call
    # if (expr.args[1] == :sin)
    #     # Pattern: sin(a*x)
    #     var_idxs = []
    #     for i in 1:length(expr.args)
    #         (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
    #         (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
    #         (expr.args[i].head == :call) && return false, expr
    #     end
    #     if length(var_idxs) == 1
    #         println("found sin term $expr")
    #         return true, expr
    #     end
    # end

    return false, expr
end

"""
    TODO: doc
"""
function resolve_cos(expr, m::PODNonlinearModel)

    function store_sin()
    end
    function lift_sin()
    end

    # @assert expr.head == :call
    # if (expr.args[1] == :cos)
    #     # Pattern: sin(a*x)
    #     var_idxs = []
    #     for i in 1:length(expr.args)
    #         (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
    #         (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
    #         (expr.args[i].head == :call) && return false, expr
    #     end
    #     if length(var_idxs) == 1
    #         println("Found cos term $expr")
    #         return true, expr
    #     end
    # end

    return false, expr
end
