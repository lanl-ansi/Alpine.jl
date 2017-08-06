"""
    expr_constr_parsing(expr, m::PODNonlinearModel)

Recognize structural constraints.
"""
function expr_constr_parsing(expr, m::PODNonlinearModel, idx::Int=0)

    isconvex = false

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.constr_patterns)
        is_strucural = eval(m.constr_patterns[i])(expr, m=m, idx=i)
        return
    end

    # Recognize built-in special structural pattern
    is_strucural = resolve_convex_constr(expr, m=m, idx=idx)

    # More patterns goes here

    return is_strucural
end

function expr_is_axn(expr, scalar=1.0, var_idxs=[]; N=nothing)

    println("inner recursive input ", expr)
    !(expr.args[1] in [:*,:^]) && return nothing, nothing         # Limited Area

    if expr.args[1] == :*
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
            elseif (expr.args[i].head == :ref)
                @assert isa(expr.args[i].args[2], Int)
                push!(var_idxs, expr.args[i])
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
                idx = expr.args[i]
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

    println("inner recursive ", scalar, " ", var_idxs)
    return scalar, var_idxs
end

"""
    A scatch for type-A convex constraint expression
"""
function resolve_convex_constr(expr; m::PODNonlinearModel=nothing, idx::Int=0)

    scalar_bin = []
    idxs_bin = []
    rhs = 0.0
    sense = expr.args[1]

    if expr.args[1] in [:(<=), :(>=)] && idx > 0
        expr_orig = :constr
        subs, rhs = expr_strip_const(expr.args[2])      # Focus on the regularized subtree (stripped with constants)
    elseif expr.args[1] == :(==)
        return false
    elseif idx == 0
        expr_orig = :obj
        subs, rhs = expr_strip_const(expr)              # Focus on the regularized subtree (stripped with constants)
    end

    @show "these are subs $subs with rhs=$rhs"
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
                    scalar, idxs = expr_is_axn(sub.args[i])
                    (scalar == nothing) && return false
                    # (mod(length(idxs),2)==0 && length(Set(idxs)) == 1) ? push!(idxs_bin, idxs) : return false
                    (length(idxs) == 2 && length(Set(idxs)) == 1) ? push!(idxs_bin, idxs[1]) : return false
                    (sub.args[1] == :-) && ((i == 2) ? push!(scalar_bin, scalar) : push!(scalar_bin, -scalar))
                    (sub.args[1] == :+) && (push!(scalar_bin, scalar))
                end
            end
        elseif sub.args[1] in [:*]
            scalar, idxs = expr_is_axn(sub, N=2)
            (scalar == nothing) && return false
            (length(idxs) == 2 && length(Set(idxs)) == 1) ? push!(idxs_bin, idxs[1]) : return false
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

    # For anything reaches this point, the return is true
    if m != nothing
        m.nonlinear_constrs[idx] = Dict(:expr_idx => idx,
                                        :nonlinear_type => :convex,
                                        :expr_orig => expr_orig,
                                        :expr_ref => expr,
                                        :convexified => false)
        m.bounding_constr_expr_mip[idx] = expr
        if expr_orig == :obj
            m.structural_obj = :convex
            m.bounding_obj_mip = Dict(:sense => nothing,
                                      :coefs => scalar_bin,
                                      :vars => idxs_bin,
                                      :rhs => -rhs,
                                      :cnt => length(idxs_bin))
        elseif expr_orig == :constr
            m.structural_constr[idx] = :convex
            m.bounding_constr_mip[idx] = Dict(:sense => sense,
                                              :coefs => scalar_bin,
                                              :vars => idxs_bin,
                                              :rhs => rhs,
                                              :cnt => length(idxs_bin))
        end
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
