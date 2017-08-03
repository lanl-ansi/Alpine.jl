# This can be a cleaner version of resolve
function expr_resolve_simple_tree(expr)
    (isa(expr, Float64) || isa(expr, Int) || isa(expr, Symbol)) && return expr
    (expr.head == :ref) && return expr

    if ((expr.args[1] == :-) && (length(expr.args) == 2))
        (isa(expr.args[2], Float64) || isa(expr.args[2], Int)) && return -expr.args[2]
        if (expr.args[2].head in [:ref, :call])
            return Expr(:call, :*, -1, expr.args[2])
        end
    end
    return expr
end

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

    function mark_convex_constr()
    end

    scalar_bin = []
    idxs_bin = []
    rhs = 0.0

    !(expr.args[1] in [:(<=), :(>=)]) && return false              # First check

    # ===================================== #
    subs, rhs = expr_strip_const(expr.args[2])   # Focus on the regularized subtree
    # ===================================== #

    @show "show ready subs $subs with rhs=$rhs"
    for sub in subs
        (isa(sub, Float64) || isa(sub, Int) || isa(sub, Symbol)) && return false
        (sub.head == :ref) && return false
        (sub.head == :call) && (length(sub.args) < 3) && return false # Second check
        if sub.args[1] in [:-, :+]
            for i in 2:length(sub.args)
                @show "OUTER $i", sub.args[i]
                if isa(sub.args[i], Float64) || isa(sub.args[i], Int)
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
        elseif sub.args[1] in [:*]
            scalar, idxs = expr_is_axn(sub, N=2)
            (scalar == nothing) && return false
            (length(idxs) == 2 && length(Set(idxs)) == 1) ? push!(idxs_bin, idxs) : return false
            push!(scalar_bin, scalar)
        else
            return false    # don't support other operators
        end

        @show rhs, scalar_bin, idxs_bin
        scalar_sign = sign(scalar_bin[1])
        if length(scalar_bin) > 1
            for i in 2:length(scalar_bin)
                !(scalar_sign == sign(scalar_bin[i])) && return false
            end
        end

        (expr.args[1] == :(<=)) && !(scalar_sign >= 0.0 && rhs >= 0.0) && return false
        (expr.args[1] == :(>=)) && !(scalar_sign <= 0.0 && rhs <= 0.0) && return false
        (expr.args[1] == :(<=)) && (scalar_sign <= 0.0) && return false
        (expr.args[1] == :(>=)) && (scalar_sign >= 0.0) && return false
    end

    return true
end

function expr_strip_const(expr, subs=[], rhs=0.0)

    exhaust_const = [!(expr.args[1] in [:+, :-]) || !(isa(expr.args[i], Float64) || isa(expr.args[i], Int)) for i in 2:length(expr.args)]
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
