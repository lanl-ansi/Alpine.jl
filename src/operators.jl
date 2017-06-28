#
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

"""
    process_expr(expr; kwargs...)

High-level warpper for processing expression with sub-tree operators
"""
function expr_(expr;kwargs...)
    # for #(all expressions: constraint & objective)
    #     while #(expr is stable)
    #         # podding expression
    #     return
    # end
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
            # Evaluate the pure number sub-tree for partial flattening
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
    for i in 1:length(m.expression_patterns)
        eval(m.expression_patterns[i])(expr, m) && return expr
    end

    # Recognization of built-in structural pattern
    resolve_bilinear(expr,m) && return expr
    resolve_monomial(expr,m) && return expr
    resolve_multilinear(expr,m) && return expr
    resolve_sin(expr,m) && return expr
    resolve_cos(expr,m) && return expr

    return expr # if no structure is detected, simply return the original tree
end

function resolve_bilinear(expr, m::PODNonlinearModel)

    @assert expr.head == :call

    function store_bilinear()
        y_idx = length(keys(nonlinear_info)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, yidx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, var_idxs[1], var_idxs[2]))
        m.nonlinear_info[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                    :ref => term_key,
                                    :lifted_constr_ref => lifted_constr_ref,
                                    :monomial_status => false,    # Later to be deprecated
                                    :cos_status => false,                   # Later to be deprecated
                                    :sin_status => false,                   # Later to be deprecated
                                    :nonlinear_type => :bilinear)
    end

    function lift_bilinear()
        expr = Expr(:*, m.nonlinear_info[term_key][:lifted_var_ref], scalar)
    end

    if (expr.args[1] == :*)  # confirm head (:*)
        # ----- Pattern A : coefficients * x * y  ------ #
        # Collect children information for checking
        scalar = 1.0
        var_idxs = []
        for i in 2:length(expr.args)
            if (isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        # Confirm patter A
        if length(var_idxs) == 2
            println("Found bilinear term $expr")
            # ----- Pattern A : Lifting ----- #
            term_key = [Expr(:ref, :x, idx) for idx in var_idxs]
            # Storage
            if (term in keys(m.nonlinear_info) || reverse(term) in keys(m.nonlinear_info))
                lift_bilinear()
            else
                store_bilinear()
                lift_bilinear()
            end
            return true
        else
            (m.log_level > 99) && println("Spitting subtree $expr for no bilinear pattern")
        end #---------- End of Pattern A ---------#
    end

    return false
end

function resolve_multilinear(expr, m::PODNonlinearModel)

    function store_multilinear()
        y_idx = length(keys(nonlinear_info)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, yidx)
        constr_block = "x[$(y_idx)]=="
        for idx in vars_idxs
            constr_block = string(constr_block, "x[$(idx)]")
            (idx < length(var_idxs)) && constr_block = string(constr_block, "*")
        end
        lifted_constr_ref = parse(constr_block)
        m.nonlinear_info[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                            :ref => term_key,
                                            :lifted_constr_ref => lifted_constr_ref,
                                            :monomial_status => false,    # Later to be deprecated
                                            :cos_status => false,                   # Later to be deprecated
                                            :sin_status => false,                   # Later to be deprecated
                                            :nonlinear_type => :multilinear)
    end

    function lift_multilinear()
        expr = Expr(:*, m.nonlinear_info[term_key][:lifted_var_ref], scalar)
    end

    @assert expr.head == :call
    if (expr.args[1] == :*)     # confirm head (:*)
        # Pattern: coefficients * x * y * z ...
        var_idxs = Set()
        scalar = 1.0
        for i in 1:length(expr.args)
            if (isa(expr.args[i], Float64) || isa(expr.args[i], Int))
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) > 2
            println("Found multilinear term $expr")
            term_key = Set()
            [push!(term_key, Expr(:ref, :x, idx)) for idx in var_idxs]
            if (term in keys(m.nonlinear_info) || reverse(term) in keys(m.nonlinear_info))
                lift_bilinear()
            else
                store_bilinear()
                lift_bilinear()
            end
            return true
        end
    end

    return false
end

function resolve_monomial(expr, m::PODNonlinearModel)

    function store_monomial()
        y_idx = length(keys(nonlinear_info)) + 1   # y is lifted var
        lifted_var_ref = Expr(:ref, :x, yidx)
        lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, vars_idxs[1], vars_idxs[1]))
        m.nonlinear_info[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                    :ref => term_key,
                                    :lifted_constr_ref => lifted_constr_ref,
                                    :monomial_status => true,               # Later to be deprecated
                                    :cos_status => false,                   # Later to be deprecated
                                    :sin_status => false,                   # Later to be deprecated
                                    :nonlinear_type => :monomial)
    end

    function lift_monomial()
        expr = m.nonlinear_info[term_key][:lifted_var_ref] #TODO: check if this will actually be replaced
    end

    @assert expr.head == :call
    if (expr.args[1] == :^) & length(expr.args) == 2
        # Pattern: (x)^(2)
        var_idxs = []
        power_scalar = 0
        for i in 2:length(expr.args)
            if (isa(expr.args[i], Float64) || isa(expr.args[i], Int))
                power_scalar += expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) == 1 && power_scalar == 2
            println("Found monomial term $expr")
            # ----- Pattern A : Lifting ----- #
            term_key = Set()
            for idx in var_idxs
                push!(term_key, Expr(:ref, :x, idx))
            end
            if (term in keys(m.nonlinear_info)
                lift_bilinear()
            else
                store_bilinear()
                lift_bilinear()
            end
            return true
        else
            (m.log_level > 99) && println("Spitting subtree $expr for no monomial pattern")
        end
    end

    return false
end

function resolve_sin(expr, m::PODNonlinearModel)

    function store_sin() end
    function lift_sin() end
    function recognize_sin() end

    @assert expr.head == :call
    if (expr.args[1] == :sin)
        # Pattern: sin(a*x)
        var_idxs = []
        for i in 1:length(expr.args)
            (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) == 1
            println("Found sin term $expr")
            return true
        end
    end

    return false
end

function resolve_cos(expr, m::PODNonlinearModel)

    function store_sin() end
    function lift_sin() end
    function recognize_sin() end

    @assert expr.head == :call
    if (expr.args[1] == :cos)
        # Pattern: sin(a*x)
        var_idxs = []
        for i in 1:length(expr.args)
            (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) == 1
            println("Found cos term $expr")
            return true
        end
    end

    return false
end
