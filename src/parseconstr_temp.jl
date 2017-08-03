using JuMP
using MathProgBase

function expr_resolve_simple_tree(expr)
    (isa(expr, Float64) || isa(expr, Int) || isa(expr, Symbol)) && return expr
    (expr.head == :ref) && return expr

    if ((expr.args[1] == :-) && (length(expr.args) == 2))
        (isa(expr.args[2], Float64) || isa(expr.args[2], Int)) && return -expr.args[2]
    end
    return expr
end


"""
    Recursivly check if a sub-tree(:call) is with a structure of
        (a * x^2) or (a * x * x)
    where a is a coefficient, and x is a variable reference (:ref)
"""
function expr_is_axn(expr, scalar=1.0, var_idxs=[]; N=nothing)
    @show "inner input ", expr

    !(expr.args[1] in [:*,:^]) && return nothing, nothing         # Limited Area

    if expr.args[1] == :*
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
            elseif (expr.args[i].head == :ref)
                @assert isa(expr.args[i].args[2], Int)
                push!(var_idxs, expr.args[i].args[2])
            elseif (expr.args[i].head == :call)
                scalar, var_idxs = expr_is_axn(expr.args[i], scalar, var_idxs)
            end
        end
    else
        power = 0
        idx = 0
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                power = expr.args[i]
                continue
            elseif (expr.args[i].head == :ref)
                idx = expr.args[i].args[2]
                continue
            elseif (expr.args[i].head == :call)
                return nothing, nothing
            end
        end
        (power <= 0) && return nothing, nothing       # Unlikely but needs check
        for i in 1:power
            push!(var_idxs, idx)
        end
    end

    !(N == nothing) && !(length(var_idxs) == N) && return nothing, nothing

    @show "inner", scalar, var_idxs
    return scalar, var_idxs
end

function resolve_convex_constr(expr)

    scalar_bin = []
    idxs_bin = []
    rhs = 0.0

    !(expr.args[1] in [:(<=), :(>=)]) && return false              # First check


    # ===================================== #
    subs, rhs = expr_strip_const(expr)   # Focus on the regularized subtree
    # ===================================== #

    @show "show ready subs $subs"
    for sub in subs

        (isa(sub, Float64) || isa(sub, Int) || isa(sub, Symbol)) && return false
        (sub.head == :ref) && return false
        (sub.head == :call) && (length(sub.args) < 3) && return false # Second check

        !(sub.args[1] in [:-, :+]) && return false              # Must be a +/- tree

        for i in 2:length(sub.args)
            @show "OUTER ", sub.args[i]
            if isa(sub.args[i], Float64) || isa(expr.args[i], Int)
                (sub.args[1] == [:-]) && ((i == 1) ? rhs += sub.args[i] : rhs -= sub.args[i])
                (sub.args[1] == [:+]) && (rhs += sub.args[i])
            elseif (sub.args[i].head == :ref)
                return false
            elseif (sub.args[i].head == :call)
                scalar, idxs = expr_is_axn(sub.args[i], N=2)
                (scalar == nothing) && return false
                (length(idxs) == 2 && length(Set(idxs)) == 1) ? push!(idxs_bin, idxs) : return false
                (sub.args[1] == :-) && ((i == 2) ? push!(scalar_bin, scalar) : push!(scalar_bin, -scalar))
                (sub.args[1] == :+) && (push!(scalar_bin, scalar))
            end
        end

        @show rhs, scalar_bin, idxs_bin
        scalar_sign = sign(scalar_bin[1])
        if length(scalar_bin) > 1
            for i in 2:length(scalar_bin)
                !(scalar_sign == sign(scalar_bin[i])) && return false
            end
        end

        (expr.args[1] == :(<=)) && (scalar_sign <= 0.0 && rhs >= 0.0) && return false
        (expr.args[1] == :(>=)) && (scalar_sign >= 0.0 && rhs <= 0.0) && return false
    end

    return true
end

function expr_strip_const(expr, subs=[], rhs=0.0)

    exhaust_const = [!(isa(expr.args[i], Float64) || isa(expr.args[i], Int)) for i in 2:length(expr.args)]
    if prod(exhaust_const)
        push!(subs, expr)
        return subs, rhs
    end

    for i in 2:length(expr.args)
        if (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol))
            (expr.args[1] == :+) && (rhs -= expr.args[i])
            (expr.args[1] == :-) && ((i == 2) ? rhs -= expr.args[i] : rhs += expr.args[i])
        elseif expr.args[i].head == :ref
            continue
        elseif expr.args[i].head == :call
            subs, rhs = expr_strip_const(expr.args[i], subs, rhs)
        end
    end

    return subs, rhs
end

function preprocess_expression(expr)

    for i in 1:length(expr.args)
        expr.args[i] = expr_resolve_simple_tree(expr.args[i])
        if (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol))
            continue
        elseif (expr.args[i].head == :ref)
            continue
        elseif (expr.args[i].head == :call)
            preprocess_expression(expr.args[i])
        end
    end

    return
end

m = Model()

@variable(m, x[1:5])
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)       # Pass
@NLconstraint(m, (3*x[1]*x[1] + 4*x[2]*x[2]) <= 25)     # Pass
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] - 25 <= 0)   # NOT Pass
@NLconstraint(m, -3*x[1]*x[1] -4*x[2]*x[2] >= -25)      # Pass
@NLconstraint(m, 3*x[1]*x[1] + 5x[2]*x[2] <= 25)        # Pass
@NLconstraint(m, 4*x[1]^2 + 5x[2]^2 <= 25)              # Pass
@NLconstraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)   # Unsupported
@NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[1] <= 25)       # Pass -> false

test_range = 8

d = JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:ExprGraph])

exs = []
exs_convex = []
for ex_i in 1:test_range
    push!(exs, MathProgBase.constr_expr(d, ex_i))
    preprocess_expression(exs[end])
    push!(exs_convex, resolve_convex_constr(exs[end]))
    println("-------")
end

println("Convexity | Expression")
for ex_i in 1:test_range
    if exs_convex[ex_i]
        println("$(exs_convex[ex_i])  | $(exs[ex_i])")
    else
        println("$(exs_convex[ex_i]) | $(exs[ex_i])")
    end
end
