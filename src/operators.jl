"""
    process_expr(expr; kwargs...)

High-level warpper for processing expression with sub-tree operators
"""
function process_expr(m::PODNonlinearModel)

    expr_initialization(m)      # S0 : initialize the space for parsing and analyzing
    expr_preprocess(m)          # S1 : pre-process the negative sign in expressions
    expr_parsing(m)             # S2 : parsing the expressions for nonlinear information
    expr_conversion(m)          # S3 : convert lifted(linear) expressions into affine function
    expr_finalized(m)           # S4 : finalize process by extracting some measurements

    return
end

"""
	STEP 1: initialize the expression/ space
"""
function expr_initialization(m::PODNonlinearModel)

	# 0 : deepcopy data into mip lifted expr place holders
	m.bounding_obj_expr_mip = deepcopy(m.obj_expr_orig)
	m.bounding_obj_mip = Dict()

	for i in 1:m.num_constr_orig
		push!(m.bounding_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
		push!(m.bounding_constr_mip, Dict())
	end

	return
end

"""
	STEP 2: preprocess expression for trivial sub-trees and nasty pieces for easier later process
"""
function expr_preprocess(m::PODNonlinearModel)

	expr_resolve_const(m.bounding_obj_expr_mip)
	expr_resolve_sign(m.bounding_obj_expr_mip)
	expr_flatten(m.bounding_obj_expr_mip)
	for i in 1:m.num_constr_orig
		expr_resolve_const(m.bounding_constr_expr_mip[i])
		expr_resolve_sign(m.bounding_constr_expr_mip[i])
		expr_flatten(m.bounding_constr_expr_mip[i].args[2])
	end

	return
end

"""
	STEP 3: parse expression for patterns on either the generic level or term level
"""
function expr_parsing(m::PODNonlinearModel)

	is_strucural = expr_constr_parsing(m.bounding_obj_expr_mip, m)
	if !is_strucural
		m.bounding_obj_expr_mip = expr_term_parsing(m.bounding_obj_expr_mip, 0, m)
		m.structural_obj = :generic_linear
	end
	(m.log_level > 99) && println("$(m.structural_obj) | [OBJ] $(m.obj_expr_orig)")

	for i in 1:m.num_constr_orig
		is_strucural = expr_constr_parsing(m.bounding_constr_expr_mip[i], m, i)
		if !is_strucural
			m.bounding_constr_expr_mip[i] = expr_term_parsing(m.bounding_constr_expr_mip[i], i, m)
			m.structural_constr[i] = :generic_linear
		end
		(m.log_level > 99) && println("$(m.structural_constr[i]) | [CONSTR] $(m.constr_expr_orig[i])")
	end

	return
end

"""
	STEP 4: convert the parsed expressions into affine-based function that can be used for adding JuMP constraints
"""
function expr_conversion(m::PODNonlinearModel)

	if m.structural_obj == :generic_linear
		m.bounding_obj_mip = expr_linear_to_affine(m.bounding_obj_expr_mip)
		m.structural_obj = :affine
	end
	m.log_level > 99 && println("type :: ", m.structural_obj)
	m.log_level > 99 && println("lifted ::", m.bounding_obj_expr_mip)
	m.log_level > 99 && println("coeffs ::", m.bounding_obj_mip[:coefs])
	m.log_level > 99 && println("vars ::", m.bounding_obj_mip[:vars])
	m.log_level > 99 && println("sense ::", m.bounding_obj_mip[:sense])
	m.log_level > 99 && println("rhs ::", m.bounding_obj_mip[:rhs])
	m.log_level > 99 && println("----------------")


	for i in 1:m.num_constr_orig
		if m.structural_constr[i] == :generic_linear
			m.bounding_constr_mip[i] = expr_linear_to_affine(m.bounding_constr_expr_mip[i])
			m.structural_constr[i] = :affine
		end
		m.log_level > 99 && println("type :: ", m.structural_constr[i])
		m.log_level > 99 && println("lifted ::", m.bounding_constr_expr_mip[i])
		m.log_level > 99 && println("coeffs ::", m.bounding_constr_mip[i][:coefs])
		m.log_level > 99 && println("vars ::", m.bounding_constr_mip[i][:vars])
		m.log_level > 99 && println("sense ::", m.bounding_constr_mip[i][:sense])
		m.log_level > 99 && println("rhs ::", m.bounding_constr_mip[i][:rhs])
		m.log_level > 99 && println("----------------")
	end

	return
end


"""
	STEP 5: collect measurements and information as needed for handly operations in the algorithm section
"""
function expr_finalized(m::PODNonlinearModel)

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
	m.num_constr_convex = length([i for i in m.structural_constr if i == :convex])

	return m
end

"""
    expr_term_parsing(expr, m::PODNonlinearModel, level=0)

Recognize and process nonlinear terms in an expression
"""
function expr_term_parsing(expr, constr_id::Int, m::PODNonlinearModel, level=0; options...)

    cnt = 0
    for node in expr.args
        cnt += 1
        if isa(node, Float64) || isa(node, Int) || isa(node, Symbol)
            continue
        elseif node.head == :call
            expr.args[cnt] = expr_term_parsing(node, constr_id, m, level+1)
        elseif node.head == :ref
            continue
        else
            error("Type issue during expression parsing. ")
        end
    end

    return expr_resolve_term_pattern(expr, constr_id, m)
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
function expr_resolve_term_pattern(expr, constr_id::Int, m::PODNonlinearModel; kwargs...)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.term_patterns)
        skip, expr = eval(m.term_patterns[i])(expr, constr_id, m)
        skip && return expr
    end

    # Recognize all built-in structural patterns
    skip, expr = resolve_bilinar_term(expr, constr_id, m)
    skip && return expr

    skip, expr = resolve_monomial_term(expr, constr_id, m)
    skip && return expr

    skip, expr = resolve_multilinear_term(expr, constr_id, m)
    skip && return expr

    # resolve_sin_term(expr, constr_id, m) && return expr
    # resolve_cos_term(expr, constr_id, m) && return expr

    return expr # if no structure is detected, simply return the original tree
end

"""
    TODO: docstring
"""
function resolve_bilinar_term(expr, constr_id::Int, m::PODNonlinearModel)

    @assert expr.head == :call

    function store_bilinear_term()
        (m.log_level) > 99 && println("found bilinear term $expr")
        y_idx = m.num_var_orig + length(keys(m.nonlinear_terms)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[2])))
        m.nonlinear_terms[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                        :id => length(keys(m.nonlinear_terms)) + 1,
                                        :ref => term_key,
                                        :orig_vars => var_idxs,
                                        :evaluator => bi(k, vec) = prod([vec[i] for i in k[:orig_vars]]),
                                        :lifted_constr_ref => lifted_constr_ref,
                                        :constr_id => Set(),
                                        :nonlinear_type => :bilinear,
                                        :convexified => false)
    end

    function lift_bilinear_term()
        push!(m.nonlinear_terms[term_key][:constr_id], constr_id)
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
        if (length(var_idxs) == 2) && length(Set(var_idxs)) == 2
            term_key = [Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[2])]
            if (term_key in keys(m.nonlinear_terms) || reverse(term_key) in keys(m.nonlinear_terms))
                (term_key in keys(m.nonlinear_terms)) ? term_key = term_key : term_key = reverse(term_key)
                return true, lift_bilinear_term()
            else
                store_bilinear_term()
                return true, lift_bilinear_term()
            end
        end
    end

    return false, expr
end

"""
    TODO: docstring
"""
function resolve_multilinear_term(expr, constr_id::Int, m::PODNonlinearModel)

    function store_multilinear_term()
        (m.log_level > 99) && println("found multilinear term $expr")
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
                                            :orig_vars => var_idxs,
                                            :evaluator => bi(k, vec) = prod([vec[i] for i in k[:orig_vars]]),
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :constr_id => Set(),
                                            :nonlinear_type => :multilinear,
                                            :convexified => false)
    end

    function lift_multilinear_term()
        push!(m.nonlinear_terms[term_key][:constr_id], constr_id)
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
            term_key = []
            for idx in var_idxs
                push!(term_key, Expr(:ref, :x, idx))
            end
            if term_key in keys(m.nonlinear_terms)
                return true, lift_multilinear_term()
            else
                store_multilinear_term()
                return true, lift_multilinear_term()
            end
        end
    elseif (expr.args[1] == :^) && length(expr.args) == 3
        # Pattern: (x)^(>2)
        var_idxs = []
        power_scalar = 0
        scalar = 1.0
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                power_scalar += expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false, expr
        end
        if length(var_idxs) == 1 && power_scalar > 2.0
            term_key = []
            for i in 1:power_scalar
                push!(term_key, Expr(:ref, :x, var_idxs[1]))
            end
            if term_key in keys(m.nonlinear_terms)
                return true, lift_multilinear_term()
            else
                store_multilinear_term()
                return true, lift_multilinear_term()
            end
        end
    end

    return false, expr
end

"""
    TODO: docstring
"""
function resolve_monomial_term(expr, constr_id::Int, m::PODNonlinearModel)

    function store_monomial_term()
        (m.log_level > 99) && println("found monomial term $expr")
        y_idx = m.num_var_orig + length(keys(m.nonlinear_terms)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[1])))
        m.nonlinear_terms[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :id => length(keys(m.nonlinear_terms)) + 1,
                                            :ref => term_key,
                                            :orig_vars => var_idxs,
                                            :evaluator => bi(k, vec) = prod([vec[i] for i in k[:orig_vars]]),
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :constr_id => Set(),
                                            :nonlinear_type => :monomial,
                                            :convexified => false)
    end

    function lift_monomial_term()
        push!(m.nonlinear_terms[term_key][:constr_id], constr_id)
        if scalar == 1.0
            return m.nonlinear_terms[term_key][:lifted_var_ref]
        else
            return Expr(:call, :*, m.nonlinear_terms[term_key][:lifted_var_ref], scalar)
        end
    end

    # Type 1 monomial
    if (expr.args[1] == :^) && length(expr.args) == 3
        # Pattern: (x)^(2)
        var_idxs = []
        power_scalar = 0
        scalar = 1.0
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
            term_key = []
            for i in 1:2
                push!(term_key, Expr(:ref, :x, var_idxs[1]))
            end

            if term_key in keys(m.nonlinear_terms)
                return true, lift_monomial_term()
            else
                store_monomial_term()
                return true, lift_monomial_term()
            end
        end
    end

    # Type 2 monomial term : x * x
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
        if (length(var_idxs) == 2) && (length(Set(var_idxs)) == 1)
            term_key = [Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[2])]
            if term_key in keys(m.nonlinear_terms)
                term_key in keys(m.nonlinear_terms)
                return true, lift_monomial_term()
            else
                store_monomial_term()
                return true, lift_monomial_term()
            end
        end
    end

    return false, expr
end

"""
    TODO: doc
"""
function resolve_sin_terme(expr, constr_id::Int, m::PODNonlinearModel)

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
function resolve_cos_term(expr, constr_id::Int, m::PODNonlinearModel)

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

"""
    A scatch for type-A convex constraint expression
"""
function resolve_convex_constr(expr, m::PODNonlinearModel=nothing, idx::Int=0, scalar_bin=[], idxs_bin=[], power_bin=[], rhs=0.0)

    if expr.args[1] in [:(<=), :(>=)] && idx > 0
        expr_orig = :constr
        sense = expr.args[1]
        rhs = expr.args[3]
        subs, rhs_strip = expr_strip_const(expr.args[2])      # Focus on the regularized subtree (stripped with constants)
        rhs += rhs_strip                                      # TODO: check if sign is flipped
    elseif expr.args[1] == :(==)
        return false
    elseif idx == 0
        expr_orig = :obj
        subs, rhs = expr_strip_const(expr)              # Focus on the regularized subtree (stripped with constants)
    end

    for sub in subs
        (isa(sub, Float64) || isa(sub, Int) || isa(sub, Symbol)) && return false
        !(sub.args[1] in [:+,:-,:*]) && return false
        (sub.head == :ref) && return false
        (sub.head == :call) && (length(sub.args) < 3) && return false
        if sub.args[1] in [:-, :+]
            for i in 2:length(sub.args)
                if isa(sub.args[i], Float64) || isa(sub.args[i], Int)
                    (sub.args[1] == [:-]) && ((i == 1) ? rhs += sub.args[i] : rhs -= sub.args[i])
                    (sub.args[1] == [:+]) && (rhs += sub.args[i])
                elseif (sub.args[i].head == :ref)
                    return false
                elseif (sub.args[i].head == :call)
                    scalar, idxs, powers = expr_is_axn(sub.args[i])
                    (scalar == nothing) && return false
                    if length(Set(idxs)) == 1
                        push!(idxs_bin, idxs[1])
                        push!(power_bin, sum(powers))
                    else
                        return false
                    end
                    (sub.args[1] == :-) && ((i == 2) ? push!(scalar_bin, scalar) : push!(scalar_bin, -scalar))
                    (sub.args[1] == :+) && (push!(scalar_bin, scalar))
                end
            end
        elseif sub.args[1] in [:*]
            scalar, idxs, powers = expr_is_axn(sub)
            (scalar == nothing) && return false
            if length(Set(idxs)) == 1
                push!(idxs_bin, idxs[1])
                push!(power_bin, sum(powers))
            else
                return false
            end
            push!(scalar_bin, scalar)
        else
            return false    # don't support other operators
        end

        # @show rhs, scalar_bin, idxs_bin
        scalar_sign = sign(scalar_bin[1])
        if length(scalar_bin) > 1
            for i in 2:length(scalar_bin)
                !(scalar_sign == sign(scalar_bin[i])) && return false
            end
        end

        !(length(Set(power_bin)) == 1) && return false
        (expr.args[1] == :(<=)) && !(scalar_sign >= 0.0 && rhs >= 0.0) && return false
        (expr.args[1] == :(>=)) && !(scalar_sign <= 0.0 && rhs <= 0.0) && return false
        (expr.args[1] == :(<=)) && (scalar_sign <= 0.0) && return false
        (expr.args[1] == :(>=)) && (scalar_sign >= 0.0) && return false
    end

    if expr_orig == :constr
        # For anything reaches this point, we've detected most common attributes of a convex expression
        # Now, we should use the differences to indicate different types of convex expression
        convex_type = :Unknown

        # Convex constraint Type-A
        # sum_i (c_i*x_i^p) <= K (K > 0, p is an even positive integer and x_i \in R)   => Convex
        power_check = [((i > 1.0) && mod(i, 2) == 0) for i in power_bin]
        (convex_type == :Unknown) && prod(power_check) && (convex_type = :convexA)

        # Type-B: Convex constraint
        # sum_i (c_i*x_i^p) <= K (K > 0, p >= 1 and x_i >= 0)                           => Convex
        lb_check = [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        power_check = [(i > 1.0) for i in power_bin]
        (convex_type == :Unknown) && prod(lb_check) && prod(power_check) && (convex_type = :convexB)

        # Convex constraint Type-C
        # sum_i (c_i*x_i^p) >= K (K > 0, 0 < p < 1 and x_i >= 0)                        => Convex
        power_check = [((i<1) && (i>0.0)) for i in power_bin]
        lb_check = [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        (convex_type == :Unknown) && prod(power_check) && prod(power_check) && (convex_type = :convexC)

        # Convex constraint Type-D
        # sum_i (c_i*x_i^p) <= K (K > 0, p <= 0 and x_i > 0)                            => Convex
        power_check = [i <= 0.0 for i in power_bin]
        lb_check =  [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        (convex_type == :Unknown) && prod(power_check) && prod(power_check) && (convex_type = :convexD)

        (convex_type == :Unknown) && return false

        # Recording the nonlinear info
        m.nonlinear_constrs[idx] = Dict(:expr_idx => idx,
                                        :convex_type => convex_type,
                                        :nonlinear_type => :convex,
                                        :expr_orig => :constraints,
                                        :expr_ref => deepcopy(expr),
                                        :convexified => false)
        # Recording the un-changed expression
        m.bounding_constr_expr_mip[idx] = expr

        # Recording structural identifier
        m.structural_constr[idx] = :convex

        # Recording the function for adding constraints
        m.bounding_constr_mip[idx] = Dict(:sense => sense,
                                          :coefs => scalar_bin,
                                          :vars => idxs_bin,
                                          :rhs => rhs,
                                          :cnt => length(idxs_bin),
                                          :powers => power_bin)
        return true

    elseif expr_orig == :obj

        convex_type = :Unknown

        # Follows the same mapping to convex constraints

        # Type-A: Convex objective function
        # sum_i (c_i*x_i^p) (p is an even positive integer and x_i \in R)   => Convex
        power_check = [((i > 1.0) && mod(i, 2) == 0) for i in power_bin]
        (convex_type == :Unknown) && prod(power_check) && (convex_type = :convexA)

        # Type-B: Convex objective function
        # sum_i (c_i*x_i^p) (p >= 1 and x_i >= 0)                           => Convex
        lb_check = [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        power_check = [i > 1.0 for i in power_bin]
        (convex_type == :Unknown) && prod(lb_check) && prod(power_check) && (convex_type = :convexB)

        # Type-D: Convex objective function
        # sum_i (c_i*x_i^p) (p <= 0 and x_i > 0) => Convex
        power_check = [i <= 0.0 for i in power_bin]
        lb_check =  [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        (convex_type == :Unknown) && prod(power_check) && prod(power_check) && (convex_type = :convexD)

        convex_type == :Unknown && return false

        # Recording the nonlinear info
        m.nonlinear_constrs[idx] = Dict(:expr_orig => :objective,
                                        :expr_idx => idx,
                                        :convex_type => convex_type,
                                        :nonlinear_type => :convex,
                                        :expr_ref => deepcopy(expr),
                                        :convexified => false)

        # Recording the un-changed expression
        m.bounding_obj_expr_mip = expr

        # Recording structural identifier
        m.structural_obj = :convex

        # Record a function that can be used to add objective function
        m.bounding_obj_mip = Dict(:sense => nothing,
                                  :coefs => scalar_bin,
                                  :vars => idxs_bin,
                                  :rhs => -rhs,
                                  :cnt => length(idxs_bin),
                                  :powers => power_bin)


    end

    return true
end
