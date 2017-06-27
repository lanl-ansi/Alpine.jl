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
function expr_podding(expr, m::PODNonlinearModel, level=0;options...)

    cnt = 0
    for node in expr.args
        cnt += 1
        if isa(node, Float64) || isa(node, Int) || isa(node, Symbol)
            continue
        elseif node.head == :call
            # if isa(eval(node), Float64) || isa(eval(node), Int)
            #     expr.args[cnt] = eval(node)
            #     continue
            # end
            expr.args[cnt] = expr_podding(node, m, level+1)
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
    if (expr.args[1] == :*)  # confirm head (:*)
        # Pattern: coefficients * x * y
        var_idxs = []
        for i in 1:length(expr.args)
            (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) == 2
            println("Found bilinear term $expr")
            return true
        end
    end

    return false
end

function resolve_multilinear(expr, m::PODNonlinearModel)

    @assert expr.head == :call
    if (expr.args[1] == :*)     # confirm head (:*)
        # Pattern: coefficients * x * y * z ...
        var_idxs = []
        for i in 1:length(expr.args)
            (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) > 2
            println("Found multilinear term $expr")
            return true
        end
    end

    return false
end

function resolve_monomial(expr, m::PODNonlinearModel)

    @assert expr.head == :call
    if (expr.args[1] == :^) & length(expr.args) == 2
        # Pattern: (x)^(2)
        var_idxs = []
        for i in 1:length(expr.args)
            (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :call) && return false
        end
        if length(var_idxs) == 1
            println("Found monomial term $expr")
            return true
        end
    end

    return false
end

function resolve_sin(expr, m::PODNonlinearModel)

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
