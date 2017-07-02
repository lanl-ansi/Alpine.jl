# =====------ Developer's Note ------===== #
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
#         / | .. \
#       x   z     y
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
# =====------ Developer's Note ------===== #

"""
    process_expr(expr; kwargs...)

High-level warpper for processing expression with sub-tree operators
"""
function expr_batch_proces(m::PODNonlinearModel;kwargs...)

	# 0 : migrate data into mip place holders
	m.lifted_obj_expr_mip = deepcopy(m.obj_expr_orig)
	for i in 1:m.num_constr_orig
        push!(m.lifted_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
    end

    # 1 : pre-process for negative
    expr_resolve_sign(m.lifted_obj_expr_mip)
    expr_flatten(m.lifted_obj_expr_mip)
    for i in 1:m.num_constr_orig
        expr_resolve_sign(m.lifted_constr_expr_mip[i])
        expr_flatten(m.lifted_constr_expr_mip[i])
    end

    # 1 : most important
    m.lifted_obj_expr_mip = expr_parsing(m.lifted_obj_expr_mip, m)
    for i in 1:m.num_constr_orig
        m.lifted_constr_expr_mip[i] = expr_parsing(m.lifted_constr_expr_mip[i], m)
        m.log_level > 99 && println(m.lifted_constr_expr_mip[end])
    end

    # 2: extract side information
    for i in keys(m.nonlinear_info)
        for var in i
            @assert isa(var.args[2], Int)
            if !(var.args[2] in m.all_nonlinear_vars)
                push!(m.all_nonlinear_vars, var.args[2])
            end
        end
    end
    m.all_nonlinear_vars = sort(m.all_nonlinear_vars)
    m.num_var_lifted_mip = length(m.nonlinear_info)

    return m
end

"""
    expr_podding(expr, m::PODNonlinearModel, level=0)
"""
function expr_parsing(expr, m::PODNonlinearModel, level=0;options...)

    cnt = 0
    for node in expr.args
        cnt += 1
        if isa(node, Float64) || isa(node, Int) || isa(node, Symbol)
            continue
        elseif node.head == :call
            # TODO: Evaluate the pure number sub-tree for partial flattening
            # if isa(eval(node), Float64) || isa(eval(node), Int)
            #     expr.args[cnt] = eval(node)
            #     continue
            # end
            expr.args[cnt] = expr_parsing(node, m, level+1)
        elseif node.head == :ref
            continue
        else
            error("Type issue during expression podding. ")
        end
    end

    return expr_resolve_pattern(expr, m)
end

"""
    expr_resolve_pattern(expr, m::PODNonlinearModel)

This function recognize, store, and replace a sub-tree `expr` with available
user-defined/built-in structures patterns. The procedure is to inqury user-defined
structure first and perform any user-defined operations to convexify the problem.
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
function expr_resolve_pattern(expr, m::PODNonlinearModel; kwargs...)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.expr_patterns)
        skip, expr = eval(m.expr_patterns[i])(expr, m)
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
    TODO: doc
"""
function resolve_bilinear(expr, m::PODNonlinearModel)

    @assert expr.head == :call

    function store_bilinear()
        y_idx = m.num_var_orig + length(keys(m.nonlinear_info)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, var_idxs[1], var_idxs[2]))
        m.nonlinear_info[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :id => y_idx,
                                            :ref => term_key,
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :monomial_status => false,    # Later to be deprecated
                                            :nonlinear_type => :bilinear,
                                            :convexified => false)
    end

    function lift_bilinear()
        if scalar == 1
            return m.nonlinear_info[term_key][:lifted_var_ref]
        else
            return Expr(:call, :*, m.nonlinear_info[term_key][:lifted_var_ref], scalar)
        end
    end

    if (expr.args[1] == :*)  # confirm head (:*)
        # ----- Pattern A : coefficients * x * y  ------ #
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
        # Confirm patter A
        if length(var_idxs) == 2
            (m.log_level) > 99 && println("found bilinear term $expr")
            term_key = [Expr(:ref, :x, idx) for idx in var_idxs]
            if (term_key in keys(m.nonlinear_info) || reverse(term_key) in keys(m.nonlinear_info))
                return true, lift_bilinear()
            else
                store_bilinear()
                return true, lift_bilinear()
            end
        end #---------- End of Pattern A ---------#
    end

    return false, expr
end

"""
    TODO: doc
"""
function resolve_multilinear(expr, m::PODNonlinearModel)

    function store_multilinear()
        y_idx = m.num_var_orig + length(keys(m.nonlinear_info)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        constr_block = "x[$(y_idx)]=="
        for j in 1:length(var_idxs)
            constr_block = string(constr_block, "x[$(var_idxs[j])]")
            if j < length(var_idxs)
                constr_block=string(constr_block, "*")
            end
        end
        lifted_constr_ref = parse(constr_block)
        m.nonlinear_info[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :id => y_idx,
                                            :ref => term_key,
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :monomial_status => false,
                                            :nonlinear_type => :multilinear,
                                            :convexified => false)
    end

    function lift_multilinear()
        if scalar == 1
            return m.nonlinear_info[term_key][:lifted_var_ref]
        else
            return Expr(:call, :*, m.nonlinear_info[term_key][:lifted_var_ref], scalar)
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
            if term_key in keys(m.nonlinear_info)
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
    TODO: doc
"""
function resolve_monomial(expr, m::PODNonlinearModel)

    function store_monomial()
        y_idx = m.num_var_orig + length(keys(m.nonlinear_info)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, y_idx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, var_idxs[1], var_idxs[1]))
        m.nonlinear_info[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :id => y_idx,
                                            :ref => term_key,
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :nonlinear_type => :monomial,
                                            :convexified => false)
    end

    function lift_monomial()
        new_expr = m.nonlinear_info[term_key][:lifted_var_ref] #TODO: check if this will actually be replaced
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

            if term_key in keys(m.nonlinear_info)
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
function resolve_sin(expr, m::PODNonlinearModel)

    function store_sin() end
    function lift_sin() end

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

    function store_sin() end
    function lift_sin() end

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
