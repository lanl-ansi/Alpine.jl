"""
    expr_constr_parsing(expr, m::PODNonlinearModel)

Recognize structural constraints.
"""
function expr_constr_parsing(expr, m::PODNonlinearModel, idx::Int=0)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.constr_patterns)
        is_strucural = eval(m.constr_patterns[i])(expr, m, idx)
        return
    end

    # Recognize built-in special structural pattern
    if m.recognize_convex
        is_convex = resolve_convex_constr(expr, m, idx)
        is_convex && return true
    end

    # More patterns goes here

    return false
end

function expr_is_axn(expr, scalar=1.0, var_idxs=[], power=[]; N=nothing)

    # println("inner recursive input ", expr)
    !(expr.args[1] in [:*,:^]) && return nothing, nothing, nothing         # Limited Area

    if expr.args[1] == :*
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
            elseif (expr.args[i].head == :ref)
                @assert isa(expr.args[i].args[2], Int)
                push!(var_idxs, expr.args[i])
                push!(power, 1)
            elseif (expr.args[i].head == :call)
                scalar, var_idxs, power = expr_is_axn(expr.args[i], scalar, var_idxs, power)
            end
        end
    elseif expr.args[1] == :^
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                push!(power, expr.args[i])
                continue
            elseif (expr.args[i].head == :ref)
                push!(var_idxs, idx)
                continue
            elseif (expr.args[i].head == :call)
                return nothing, nothing, nothing
            end
        end
    end

    # If the user wants a specific N
    !(N == nothing) && !(length(var_idxs) == N) && return nothing, nothing, nothing

    # println("inner recursive ", scalar, " ", var_idxs)
    @assert length(var_idxs) == length(power)
    return scalar, var_idxs, power
end

"""
    A scatch for type-A convex constraint expression
"""
function resolve_convex_constr(expr, m::PODNonlinearModel=nothing, idx::Int=0)

    scalar_bin = []
    idxs_bin = []
    power_bin = []
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

        # Type-A: Convex constraint
        # sum_i (c_i*x_i^p) <= K (K > 0, p >= 1 and x_i >= 0)                           => Convex
        lb_check = [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        power_check = [(i >= 1.0) for i in power_bin]
        prod(lb_check) && prod(power_check) && (convex_type = :convexA)

        # Convex constraint Type-B
        # sum_i (c_i*x_i^p) <= K (K > 0, p is an even positive integer and x_i \in R)   => Convex
        power_check = [((i > 1.0) && mod(i, 2) == 0) for i in power_bin]
        prod(power_check) && (convex_type = :convexB)

        # Convex constraint Type-C
        # sum_i (c_i*x_i^p) >= K (K > 0, 0 < p < 1 and x_i >= 0)                        => Convex
        power_check = [((i<1) && (i>0.0)) for i in power_bin]
        lb_check = [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        prod(power_check) && prod(power_check) && (convex_type = :convexC)

        # Convex constraint Type-D
        # sum_i (c_i*x_i^p) <= K (K > 0, p <= 0 and x_i > 0)                            => Convex
        power_check = [i <= 0.0 for i in power_bin]
        lb_check =  [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        prod(power_check) && prod(power_check) && (convex_type = :convexD)

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
        # sum_i (c_i*x_i^p) (p >= 1 and x_i >= 0)                           => Convex
        lb_check = [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        power_check = [i >= 1.0 for i in power_bin]
        prod(lb_check) && prod(power_check) && (convex_type = :convexA)

        # Type-B: Convex objective function
        # sum_i (c_i*x_i^p) (p is an even positive integer and x_i \in R)   => Convex
        power_check = [((i > 1.0) && mod(i, 2) == 0) for i in power_bin]
        prod(power_check) && (convex_type = :convexB)

        # Type-D: Convex objective function
        # sum_i (c_i*x_i^p) (p <= 0 and x_i > 0) => Convex
        power_check = [i <= 0.0 for i in power_bin]
        lb_check =  [m.l_var_orig[i.args[2]] >= 0.0 for i in idxs_bin]
        prod(power_check) && prod(power_check) && (convex_type = :convexD)

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
