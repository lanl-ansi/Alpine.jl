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

    return detect_nonconvex_terms(expr, constr_id, m)
end

"""
    detect_nonconvex_terms(expr, m::PODNonlinearModel)

This function recognizes, stores, and replaces a sub-tree `expr` with available
user-defined/built-in structures patterns. The procedure is creates the required number
of lifted variables based on the patterns that it it trying to recognize.
Then, go through all built-in structures and perform operatins to convexify the problem.

Specific structure pattern information will be described formally.
"""
function detect_nonconvex_terms(expr, constr_id::Int, m::PODNonlinearModel; kwargs...)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.term_patterns)
        skip, expr = eval(m.term_patterns[i])(expr, constr_id, m)
        skip && return expr
    end

    # NOTE :: Sequence of which term to detect matters here
    #         More specific term should be detected first to reduce the size of relaxation
    #         model. For example, power-2 or bilinear term should be detected before
    #         general multilinear terms to apply better outter approxiamtion.

    # LEVEL 1 : : Recognize all built-in structural patterns
    skip, expr = detect_discretemulti_term(expr, constr_id, m)     #L1 : General case : discrete variables * continous variables
    skip && return expr

    skip, expr = detect_binprod_term(expr, constr_id, m)           #L1 : binprod = binary products
    skip && return expr

    skip, expr = detect_intprod_term(expr, constr_id, m)		   #L1 : intprod = integer products
    skip && return expr

    # LEVEL 2 : : Recognize all built-in structural patterns
    skip, expr = detect_bilinear_term(expr, constr_id, m)          #L2
    skip && return expr

    skip, expr = detect_monomial_term(expr, constr_id, m)          #L2 : must go before multilinear
    skip && return expr

    skip, expr = detect_multilinear_term(expr, constr_id, m)       #L2
    skip && return expr

    skip, expr = detect_sincos_term(expr, constr_id, m)            #L2
    skip && return expr

    return expr # if no structure is detected, simply return the original tree
end

function store_nonconvex_term(m::PODNonlinearModel, nl_key::Any, var_idxs::Any, term_type::Symbol, operator::Symbol, evaluator::Function, bd_resolver::Function, discvar_collector::Function)

    l_cnt = length(keys(m.linear_terms))
    nl_cnt = length(keys(m.nonconvex_terms))

    y_idx = m.num_var_orig + nl_cnt + l_cnt  + 1   # y is always lifted var
    lifted_var_ref = Expr(:ref, :x, y_idx)
    lifted_constr_ref = build_constr_block(y_idx, var_idxs, operator)

    m.nonconvex_terms[nl_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                    :id => nl_cnt + 1,
                                    :y_idx => y_idx,
                                    :y_type => resolve_lifted_var_type([m.var_type[k] for k in var_idxs], operator),
                                    :var_idxs => var_idxs,
                                    :ref => nl_key,
                                    :evaluator => evaluator,
                                    :lifted_constr_ref => lifted_constr_ref,
                                    :constr_id => Set(),
                                    :nonlinear_type => term_type,
                                    :convexified => false,
                                    :bound_resolver => bd_resolver,
                                    :discvar_collector => discvar_collector)

    m.term_seq[nl_cnt+l_cnt+1] = nl_key                              # Assistive information

    # push!(m.var_type, :Cont)  # TODO check if this replacement is good since additional constraints should be able to sufficiently constraint the type
    push!(m.var_type, m.nonconvex_terms[nl_key][:y_type])            # Keep track of the lifted var type
    m.loglevel > 199 && println("found lifted $(term_type) term $(lifted_constr_ref)")
    return y_idx
end

function store_linear_term(m::PODNonlinearModel, term_key::Any, expr::Any)#, bound_resolver::Function)

    l_cnt = length(keys(m.linear_terms))
    nl_cnt = length(keys(m.nonconvex_terms))

    y_idx = m.num_var_orig + nl_cnt + l_cnt + 1

    lifted_var_ref = Expr(:ref, :x, y_idx)
    lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, expr)

    m.linear_terms[term_key] = Dict(:lifted_var_ref => lifted_var_ref,
                                    :id => length(keys(m.linear_terms)) + 1,
                                    :ref => term_key,
                                    :y_idx => y_idx,
                                    :y_type => resolve_lifted_var_type([m.var_type[k[2]] for k in term_key[:coef_var]], :+),
                                    :evaluator => linear,
                                    :lifted_constr_ref => lifted_constr_ref,
                                    :constr_id => Set(),
                                    :bound_resolver => nothing)

    m.term_seq[l_cnt+nl_cnt + 1] = term_key
    push!(m.var_type, m.linear_terms[term_key][:y_type]) # Keep track of the lifted var type
    m.loglevel > 199 && println("found lifted linear term $(lifted_var_ref) = $expr")

    return y_idx
end

function lift_nonconvex_term(m::PODNonlinearModel, nl_key, constr_id::Int, scalar = 1.0)
    push!(m.nonconvex_terms[nl_key][:constr_id], constr_id)

    if scalar == 1.0
        return m.nonconvex_terms[nl_key][:lifted_var_ref]
    else
        return Expr(:call, :*, m.nonconvex_terms[nl_key][:lifted_var_ref], scalar)
    end
end

function lift_linear_term(m::PODNonlinearModel, term_key, constr_id::Int)

    push!(m.linear_terms[term_key][:constr_id], constr_id)
    return m.linear_terms[term_key][:lifted_var_ref]

    return
end

function detect_linear_term(expr, constr_id::Int, m::PODNonlinearModel)

    @assert expr.head == :call
    coef_fetch = Dict(:+ => 1.0, :- => -1.0)

    # Re-process the expression sub-tree [REQUIRED]
    expr_resolve_const(expr)
    expr_resolve_sign(expr)
    expr_flatten(expr)

    if expr.args[1] in [:+, :-]
        scalar = 0.0
        coef_var = Set()
        for i in 2:length(expr.args)  # TODO: Consider recursive operation
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                (i == 2) ? scalar=expr.args[i] : scalar+=coef_fetch[expr.args[1]]*expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue #should never happen
            if (expr.args[i].head==:ref) && isa(expr.args[i].args[2], Int)
                (i == 2) ? push!(coef_var, (1.0, expr.args[i].args[2])) : push!(coef_var, (coef_fetch[expr.args[1]],expr.args[i].args[2]))
                continue
            end
            # Specical Check
            if expr.args[i].head == :call && expr.args[i].args[1] == :* && length(expr.args[i].args) == 3
                sub_coef = [j for j in expr.args[i].args if (isa(j, Int) || isa(j, Float64))]
                sub_vars = [j.args[2] for j in expr.args[i].args if ((:head in fieldnames(j)) && j.head == :ref)]
                (isempty(sub_coef) || isempty(sub_vars)) && return false, expr
                (length(sub_coef) != 1 || length(sub_vars) != 1) && return false, expr
                (i == 2) ? push!(coef_var, (1.0*sub_coef[1], sub_vars[1])) : push!(coef_var, (coef_fetch[expr.args[1]]*sub_coef[1], sub_vars[1]))
                continue
            end
            # General Check
            if expr.args[i].head == :call && expr.args[i].args[1] in [:-,:+] && length(expr.args[i].args) == 2 # resolve -(2) or -(x) terms
                expr.args[i].args[1] == :+ ? sub_coef = [1.0] : sub_coef = [-1.0]
                sub_vars = [j.args[2] for j in expr.args[i].args if ((:head in fieldnames(j)) && j.head == :ref)]
                isempty(sub_vars) && return false, expr
                length(sub_vars) != 1 && return false, expr
                (i == 2) ? push!(coef_var, (1.0*sub_coef[1], sub_vars[1])) : push!(coef_var, (coef_fetch[expr.args[1]]*sub_coef[1], sub_vars[1]))
            else
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                down_check ? expr.args[i] = linear_lift_var : return false, expr
                push!(coef_var, (1.0, expr.args[i].args[2]))
            end
        end
        # By reaching here, it is already certain that we have found the term, always treat with :+
        term_key = Dict(:scalar=>scalar, :coef_var=>coef_var, :sign=>:+)
        term_key in keys(m.linear_terms) || store_linear_term(m, term_key, expr)
        return true, lift_linear_term(m, term_key, constr_id)
    elseif expr.args[1] in [:*] && length(expr.args) == 3
        # For terms like (3*x)*y
        scalar = 0.0
        coef_var = Set()

        sub_coef = [i for i in expr.args if (isa(i, Int) || isa(i, Float64))]
        sub_vars = [i.args[2] for i in expr.args if ((:head in fieldnames(i)) && i.head == :ref)]
        (isempty(sub_coef) || isempty(sub_vars)) && return false, expr
        (length(sub_coef) != 1 || length(sub_vars) != 1) && return false, expr
        push!(coef_var, (1.0*sub_coef[1], sub_vars[1]))

        # By reaching here, it is already certain that we have found the term, always treat with :+
        term_key = Dict(:scalar=>scalar, :coef_var=>coef_var, :sign=>:+)
        term_key in keys(m.linear_terms) || store_linear_term(m, term_key, expr)
        return true, lift_linear_term(m, term_key, constr_id)
    end

    return false, expr
end

function basic_linear_bounds(m::PODNonlinearModel, k::Any, linear_terms=nothing)

    linear_terms == nothing ? linear_terms = m.linear_terms : linear_term = linear_terms

    lifted_idx = linear_terms[k][:y_idx]
    ub = 0.0
    lb = 0.0
    for j in linear_terms[k][:ref][:coef_var]
        (j[1] > 0.0) ? ub += abs(j[1])*m.u_var_tight[j[2]] : ub -= abs(j[1])*m.l_var_tight[j[2]]
        (j[1] > 0.0) ? lb += abs(j[1])*m.l_var_tight[j[2]] : lb -= abs(j[1])*m.u_var_tight[j[2]]
    end
    lb += linear_terms[k][:ref][:scalar]
    ub += linear_terms[k][:ref][:scalar]
    if lb > m.l_var_tight[lifted_idx] + m.tol
        m.l_var_tight[lifted_idx] = lb
    end
    if ub < m.u_var_tight[lifted_idx] - m.tol
        m.u_var_tight[lifted_idx] = ub
    end

    return
end

function basic_linear_bounds(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.linear_terms[k][:y_idx]
    ub = 0.0
    lb = 0.0
    for j in m.linear_terms[k][:ref][:coef_var]
        (j[1] > 0.0) ? ub += abs(j[1])*d[j[2]][end] : ub -= abs(j[1])*d[j[2]][1]
        (j[1] > 0.0) ? lb += abs(j[1])*d[j[2]][1] : lb -= abs(j[1])*d[j[2]][end]
    end
    lb += m.linear_terms[k][:ref][:scalar]
    ub += m.linear_terms[k][:ref][:scalar]

    if lb > d[lifted_idx][1] + m.tol
        d[lifted_idx][1] = lb
    end
    if ub < d[lifted_idx][end] - m.tol
        d[lifted_idx][end] = ub
    end

    return d
end


## Evaluators ##
bpml(k,vec) = prod([vec[i] for i in k[:var_idxs]])
discretemulti(k, vec) = prod([vec[i] for i in k[:var_idxs]])
binint(l, vec) = prod([vec[i] for i in k[:var_idxs]])
intprod(k, vec) = prod([vec[i] for i in k[:var_idxs]])
binprod(k, vec) = prod([vec[i] for i in k[:var_idxs]])
binlin(k,vec) = prod([vec[i] for i in k[:var_idxs]])
intlin(k,vec) = prod([vec[i] for i in k[:var_idxs]])
bilinear(k,vec) = prod([vec[i] for i in k[:var_idxs]])
multilinear(k,vec) = prod([vec[i] for i in k[:var_idxs]])
monomial(k, vec) = vec[k[:var_idxs][1]]^2
sincos(k, vec) = eval(k[:nonlinear_type])(vec[k[:var_idxs][1]])
linear(k, vec) = sum([i[1]*vec[i[2]] for i in k[:ref][:coef_var]])

"""
    Recognize prodcuts of binary variables and multilinear products

        General Case : x1 * x2 .. * xN * y1 * y2 .. * yM
            where x are binary variables, y are continous variables

    Leads to BINLIN terms, with BINPROD, INTPROD, INTLIN if necessary
"""
function detect_discretemulti_term(expr, constr_id::Int, m::PODNonlinearModel)

    # Alwasy construct the binlin term after lifting
    @assert expr.head == :call

    if (expr.args[1] == :*)
        # Pattern: coefficients * x * y * z ..., where x, y, z are all binary variables
        var_idxs = []
        var_types = []
        scalar = 1.0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_types, m.var_type[expr.args[i].args[2]])
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                push!(var_types, m.var_type[linear_lift_var.args[2]])
                continue
            end
        end

        term_key = [Expr(:ref, :x, idx) for idx in var_idxs]

        cont_var_idxs = [idx for idx in var_idxs if m.var_type[idx] == :Cont]
        bin_var_idxs = [idx for idx in var_idxs if m.var_type[idx] == :Bin]
        int_var_idxs = [idx for idx in var_idxs if m.var_type[idx] == :Int]

        # Must carry at both discrete and continous variables
        isempty(int_var_idxs) && isempty(bin_var_idxs) && return false, expr
        isempty(cont_var_idxs) && return false, expr

        # Lift clusters of continous variables multiplication if necessary
        if length(cont_var_idxs) > 1
            ml_term_key = [Expr(:ref, :x, idx) for idx in cont_var_idxs]
            ml_term_expr = Expr(:call, :*)
            for idx in cont_var_idxs
                push!(ml_term_expr.args, Expr(:ref, :x, idx))
            end
            ml_lift_term = detect_nonconvex_terms(ml_term_expr, constr_id, m)
            ml_idx = ml_lift_term.args[2]
        else
            ml_idx = cont_var_idxs[1]
        end

        # Lift clusters of binary vars multiplication if necessary
        if length(bin_var_idxs) > 1
            bp_term_key = [Expr(:ref, :x, idx) for idx in bin_var_idxs]
            bp_term_expr = Expr(:call, :*)
            for idx in bin_var_idxs
                push!(bp_term_expr.args, Expr(:ref, :x, idx))
            end
            bp_lift_term = detect_nonconvex_terms(bp_term_expr, constr_id, m)
            bp_idx = bp_lift_term.args[2]
        else
            isempty(bin_var_idxs) ? bp_idx = -1 : bp_idx = bin_var_idxs[1]
        end

        # lift clusters of integer varis multiplication if necessary
        if length(int_var_idxs) > 1
            ip_term_key = [Expr(:ref, :x, idx) for idx in int_var_idxs]
            ip_term_expr = Expr(:call, :*)
            for idx in int_var_idxs
                push!(ip_term_expr.args, Expr(:ref, :x, idx))
            end
            ip_lift_term = detect_nonconvex_terms(ip_term_expr, constr_id, m)
            ip_idx = ip_lift_term.args[2]
        else
            isempty(int_var_idxs) ? ip_idx = -1 : ip_idx = int_var_idxs[1]
        end

        if bp_idx < 0 # intlin term if no binary variable
            intlin_key = [Expr(:ref, :x, ip_idx), Expr(:ref, :x, ml_idx)]
            intlin_idxs = [ip_idx, ml_idx;]
            intlin_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, intlin_key, intlin_idxs, :INTLIN, :*, intlin, basic_intlin_bounds, collect_intlin_discvar)
            return true, lift_nonconvex_term(m, intlin_key, constr_id, scalar)
        end

        if ip_idx < 0 # binlin term if no integer variable
            binlin_key = [Expr(:ref, :x, bp_idx), Expr(:ref, :x, ml_idx)]
            binlin_idxs = [bp_idx, ml_idx;]
            binlin_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, binlin_key, binlin_idxs, :BINLIN, :*, binlin, basic_binlin_bounds, collect_binlin_discvar)
            return true, lift_nonconvex_term(m, binlin_key, constr_id, scalar)
        end

        # Otherwise, first :intlin term then :binlin term by assumption
        intlin_key = [Expr(:ref, :x, ip_idx), Expr(:ref, :x, ml_idx)]
        intlin_expr = Expr(:call, :*)
        push!(intlin_expr.args, Expr(:ref, :x, ip_idx))
        push!(intlin_expr.args, Expr(:ref, :x, ml_idx))
        intlin_lift = detect_nonconvex_terms(intlin_expr, constr_id, m)
        intlin_idx = intlin_lift.args[2]

        binlin_key = [Expr(:ref, :x, bp_idx), Expr(:ref, :x, intlin_idx)]
        binlin_idxs = [bp_idx, intlin_idx;]
        binlin_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, binlin_key, binlin_idxs, :BINLIN, :*, binlin, basic_binlin_bounds, collect_binlin_discvar)
        return true, lift_nonconvex_term(m, binlin_key, constr_id, scalar)
    end

    return false, expr
end

function basic_binlin_bounds(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:y_idx]

    prod_idxs = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Cont]
    @assert length(prod_idxs) == 1
    lin_idx = prod_idxs[1]

    if m.l_var_tight[lin_idx] > 0.0
        m.l_var_tight[lifted_idx] = 0.0
        m.u_var_tight[lifted_idx] = m.u_var_tight[lin_idx]
    elseif m.u_var_tight[lin_idx] < 0.0
        m.u_var_tight[lifted_idx] = 0.0
        m.l_var_tight[lifted_idx] = m.l_var_tight[lin_idx]
    else
        m.u_var_tight[lifted_idx] = m.u_var_tight[lin_idx]
        m.l_var_tight[lifted_idx] = m.l_var_tight[lin_idx]
    end

    return
end

function basic_binlin_bounds(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[k][:y_idx]

    prod_idxs = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Cont]
    @assert length(prod_idxs) == 1
    lin_idx = prod_idxs[1]

    if d[lin_idx][1] > 0.0
        d[lifted_idx][1] = 0.0
        d[lifted_idx][end] = d[lin_idx][end]
    elseif d[lin_idx][end] < 0.0
        d[lifted_idx][end] = 0.0
        d[lifted_idx][1] = d[lin_idx][1]
    else
        d[lifted_idx][end] = d[lin_idx][end]
        d[lifted_idx][1] = d[lin_idx][1]
    end

    return d
end

function basic_intlin_bounds(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:y_idx]

    lins = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Cont]
    ints = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Int]
    @assert length(lins) == 1 && length(ints) == 1

    lin_idx = lins[1]
    int_idx = ints[1]

    linrange = [m.l_var_tight[lin_idx],m.u_var_tight[lin_idx];]
    intrange = [m.l_var_tight[int_idx],m.u_var_tight[int_idx];]

    crossrange = linrange * intrange'

    m.l_var_tight[lifted_idx] = max(m.l_var_tight[lifted_idx], minimum(crossrange))
    m.u_var_tight[lifted_idx] = min(m.u_var_tight[lifted_idx], maximum(crossrange))

    return
end

function basic_intlin_bounds(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[k][:y_idx]

    lins = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Cont]
    ints = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Int]
    @assert length(lins) == 1 == length(ints)

    lin_idx = lins[1]
    int_idx = ints[1]

    linrange = [d[lin_idx][1],d[lin_idx][end];]
    intrange = [d[int_idx][1],d[int_idx][end];]

    crossrange = linrange * intrange'

    d[lifted_idx][1] = max(d[lifted_idx][1], minimum(crossrange))
    d[lifted_idx][end] = min(d[lifted_idx][end], maximum(crossrange))

    return d
end

function collect_binlin_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)
    # Exact linearization exist
    return
end

function collect_intlin_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)

    for i in m.nonconvex_terms[k][:var_idxs]
        @assert isa(i, Int)
        if var_bowl == nothing
            i in m.candidate_disc_vars || push!(m.candidate_disc_vars, i)
        else
            i in var_bowl || push!(var_bowl, i)
        end
    end

    return
end

"""
	Recognize products of discrete variables with both binary and integer variables

		General Case : x1 * x2 * .. * xN * y1 * y2 * .. * yM
			where x are binary variables, y are integer variables
		Special Case : x * y1 * .. * yM | x1 * .. * xN * y | x * y

	Leads to BININT terms, with BINPROD, INTPROD if necessary
"""
function detect_binint_term(expr, constr_id::Int, m::PODNonlinearModel)

    @assert expr.head == :call

    if (expr.args[1] == :*)

        var_idxs = []
        var_types = []
        scalar = 1.0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_types, m.var_type[expr.args[i].args[2]])
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                push!(var_types, m.var_type[linear_lift_var.args[2]])
                continue
            end
        end

        cont_var_idxs = [idx for idx in var_idxs if m.var_type[idx] == :Cont]
        bin_var_idxs = [idx for idx in var_idxs if m.var_type[idx] == :Bin]
        int_var_idxs = [idx for idx in var_idxs if m.var_type[idx] == :Int]

		# Avoid the presense of continous variables
		isempty(cont_var_idxs) || return false, expr
        isempty(int_var_idxs) && isempty(bin_var_idxs) && return false, expr

        # Lift clusters of binary vars multiplication if necessary
        if length(bin_var_idxs) > 1
            bp_term_key = [Expr(:ref, :x, idx) for idx in bin_var_idxs]
            bp_term_expr = Expr(:call, :*)
            for idx in bin_var_idxs
                push!(bp_term_expr.args, Expr(:ref, :x, idx))
            end
            bp_lift_term = detect_nonconvex_terms(bp_term_expr, constr_id, m)
            bp_idx = bp_lift_term.args[2]
        else
            bp_idx = bin_var_idxs[1]
        end

        # lift clusters of integer varis multiplication if necessary
        if length(int_var_idxs) > 1
            ip_term_key = [Expr(:ref, :x, idx) for idx in int_var_idxs]
            ip_term_expr = Expr(:call, :*)
            for idx in cont_var_idxs
                push!(ip_term_expr.args, Expr(:ref, :x, idx))
            end
            ip_lift_term = detect_nonconvex_terms(ip_term_expr, constr_id, m)
        else
            ip_idx = ip_var_idxs[1]
        end

        binint_key = [Expr(:ref, :x, bp_idx), Expr(:ref, :x, ip_idx)]
        binint_idxs = [bp_idx, ip_idx;]
        binint_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, binint_key, binint_idxs, :BININT, :*, binint, basic_binint_bound, collect_binint_discvar)
        return true, lift_nonconvex_term(m, binint_key, constr_id, scalar)
    end

    return false, expr
end

function basic_binint_bound(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:y_idx]

    prod_idxs = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Int]
    @assert length(prod_idxs) == 1
    lin_idx = prod_idxs[1]

    if m.l_var_tight[lin_idx] > 0.0
        m.l_var_tight[lifted_idx] = 0.0
        m.u_var_tight[lifted_idx] = m.u_var_tight[lifted_idx]
    elseif m.u_var_tight[lin_idx] < 0.0
        m.u_var_tight[lifted_idx] = 0.0
        m.l_var_tight[lifted_idx] = m.l_var_tight[lin_idx]
    else
        m.u_var_tight[lifted_idx] = m.u_var_tight[lin_idx]
        m.l_var_tight[lifted_idx] = m.l_var_tight[lin_idx]
    end

    return
end

function basic_binint_bound(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[k][:y_idx]

    prod_idxs = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Int]
    @assert length(prod_idxs) == 1
    lin_idx = prod_idxs[1]

    if d[lin_idx][1] > 0.0
        d[lifted_idx][1] = 0.0
        d[lifted_idx][end] = d[lin_idx][end]
    elseif m.u_var_tight[lin_idx] < 0.0
        d[lifted_idx][end] = 0.0
        d[lifted_idx][1] = d[lin_idx][1]
    else
        d[lifted_idx][end] = d[lin_idx][end]
        d[lifted_idx][1] = d[lin_idx][1]
    end

    return d
end

function collect_binint_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)
    # Exact linearization exist
    return
end

"""
    Recognize products of binary variables : x1 * x2 * .. * xN
"""
function detect_intprod_term(expr, constr_id::Int, m::PODNonlinearModel)

	@assert expr.head == :call
    if (expr.args[1] == :*)
        # Pattern: coefficients * x * y * z ..., where x, y, z are all integer variables
        var_idxs = []
        scalar = 1.0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            !isempty(var_idxs) && m.var_type[var_idxs[end]] != :Int && return false, expr
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                m.var_type[linear_lift_var.args[2]] != :Int && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                continue
            end
        end

        if length(var_idxs) >= 2
            term_key = [Expr(:ref, :x, idx) for idx in var_idxs]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :INTPROD, :*, intprod, basic_intprod_bounds, collect_intprod_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end

    elseif expr.args[1] == :^ && length(expr.args) == 3
        # Pattern: (x)^(>2), where x is integer variable
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
            !isempty(var_idxs) && m.var_type[var_idxs[end]] != :Int && return false, expr
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                m.var_type[linear_lift_var.args[2]] != :Int && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                continue
            end
        end

        if length(var_idxs) == 1 && power_scalar >= 2.0 && mod(power_scalar, 1.0) == 0.0
            term_key = [Expr(:ref, :x, var_idxs[1]) for i in 1:power_scalar]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :INTPROD, :*, intprod, basic_intprod_bounds, collect_intprod_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    end

    return false, expr
end

function basic_intprod_bounds(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]

    bound = []
    for (cnt,varexpr) in enumerate(k)
        var = varexpr.args[2]
        var_bounds = [m.l_var_tight[var], m.u_var_tight[var]]
        if cnt == 1
            bound = copy(var_bounds)
        elseif cnt == 2
            bound = bound * var_bounds'
        else
            bound = diag(bound) * var_bounds'
        end
    end

    if minimum(bound) > m.l_var_tight[lifted_idx] + m.tol
        m.l_var_tight[lifted_idx] = minimum(bound)
    end
    if maximum(bound) < m.u_var_tight[lifted_idx] - m.tol
        m.u_var_tight[lifted_idx] = maximum(bound)
    end

    return
end

function basic_intprod_bounds(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]

    bound = []
    for (cnt,varexpr) in enumerate(k)
        var = varexpr.args[2]
        var_bounds = [d[var][1], d[var][end]]
        if cnt == 1
            bound = copy(var_bounds)
        elseif cnt == 2
            bound = bound * var_bounds'
        else
            bound = diag(bound) * var_bounds'
        end
    end
    if minimum(bound) > d[lifted_idx][1] + m.tol
        d[lifted_idx][1] = minimum(bound)
    end
    if maximum(bound) < d[lifted_idx][end] - m.tol
        d[lifted_idx][end] = maximum(bound)
    end

    return d
end

function collect_intprod_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)
    for var in m.nonconvex_terms[k][:var_idxs]
        @assert isa(var, Int)
        if var_bowl == nothing
            var in m.candidate_disc_vars || push!(m.candidate_disc_vars, var)
        else
            var in var_bowl || push!(var_bowl, var)
        end
    end
    return
end

"""
    Recognize products of binary variables : x1 * x2 * .. * xN
"""
function detect_binprod_term(expr, constr_id::Int, m::PODNonlinearModel)

	@assert expr.head == :call
    if (expr.args[1] == :*)
        # Pattern: coefficients * x * y * z ..., where x, y, z are all binary variables
        var_idxs = []
        scalar = 1.0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            !isempty(var_idxs) && m.var_type[var_idxs[end]] != :Bin && return false, expr
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                m.var_type[linear_lift_var.args[2]] != :Bin && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                continue
            end
        end
        if length(var_idxs) >= 2
            term_key = [Expr(:ref, :x, idx) for idx in var_idxs]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :BINPROD, :*, binprod, basic_binprod_bounds, collect_binprod_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    elseif (expr.args[1] == :^) && length(expr.args) == 3
        # Pattern: (x)^(>2), where x is binary variable
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
            !isempty(var_idxs) && m.var_type[var_idxs[end]] != :Bin && return false, expr
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                m.var_type[linear_lift_var.args[2]] != :Bin && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                continue
            end
        end
        if length(var_idxs) == 1 && power_scalar > 2.0 && mod(power_scalar, 1.0) == 0.0
            term_key = [Expr(:ref, :x, var_idxs[1]) for i in 1:power_scalar]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :BINPROD, :*, binprod, basic_binprod_bounds, collect_binprod_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    end

    return false, expr
end

function basic_binprod_bounds(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]
    m.l_var_tight[lifted_idx] = 0
    m.u_var_tight[lifted_idx] = 1

    return
end

function basic_binprod_bounds(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]
    d[lifted_idx][1] = 0
    d[lifted_idx][end] = 1

    return d
end

function collect_binprod_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)
    # Exact linearization exists
    return
end


"""
    Future MONOMIAL Cluster
    Recognize bilinear terms: x * y, where x and y are continous variables
    Recognize multilinear terms: x1 * x2 * .. * xN, where all x are continous variables
    Recognize monomial terms: x^2 or x * x, where x is continuous
"""
function detect_bilinear_term(expr, constr_id::Int, m::PODNonlinearModel)

    @assert expr.head == :call
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
            !isempty(var_idxs) && m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr    # Don't consider the discrete variable
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr  # Don't consider the discrete variable
                continue
            end
        end

        # Cofirm detection of patter A and perform store & lifting procedures
        if (length(var_idxs) == 2) && length(Set(var_idxs)) == 2
            term_key = [Expr(:ref, :x, var_idxs[1]), Expr(:ref, :x, var_idxs[2])]
            if term_key in keys(m.nonconvex_terms) || reverse(term_key) in keys(m.nonconvex_terms)
                term_key in keys(m.nonconvex_terms) ? term_key = term_key : term_key = reverse(term_key)
                return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
            else
                store_nonconvex_term(m, term_key, var_idxs, :BILINEAR, :*, bilinear, basic_monomial_bounds, collect_monomial_discvar)
                return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
            end
        end
    end

    return false, expr
end

function detect_multilinear_term(expr, constr_id::Int, m::PODNonlinearModel)

    @assert expr.head == :call
    if (expr.args[1] == :*) # Pattern: coefficients * x * y * z ...
        var_idxs = []
        scalar = 1.0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            !isempty(var_idxs) && m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr
                continue
            end
        end
        if length(var_idxs) > 2
            term_key = [Expr(:ref, :x, idx) for idx in var_idxs]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :MULTILINEAR, :*, multilinear, basic_monomial_bounds, collect_monomial_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    elseif (expr.args[1] == :^) && length(expr.args) == 3 # Pattern: (x)^(>2)
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
            !isempty(var_idxs) && m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr # Avoid fetching terms with discrete variables
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr
                continue
            end
        end
        if length(var_idxs) == 1 && power_scalar > 2.0 && mod(power_scalar, 1.0) == 0.0
            term_key = [Expr(:ref, :x, var_idxs[1]) for i in 1:power_scalar]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :MULTILINEAR, :*, multilinear, basic_monomial_bounds, collect_monomial_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    end

    return false, expr
end

function detect_monomial_term(expr, constr_id::Int, m::PODNonlinearModel)

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
            !isempty(var_idxs) && m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr # Avoid fetching terms with discrete variables
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr # Avoid fetching terms with discrete variables
                continue
            end
        end
        if length(var_idxs) == 1 && power_scalar == 2.0
            term_key = [Expr(:ref, :x, var_idxs[1]) for i in 1:2]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :MONOMIAL, :*, monomial, basic_monomial_bounds, collect_monomial_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    end

    # Type 2 monomial term : x * x
    if (expr.args[1] == :*)  # confirm head (:*)
        # ----- Pattern : coefficient * x * y  ------ #
        scalar = 1.0
        var_idxs = []
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            (isa(expr.args[i], Symbol)) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            !isempty(var_idxs) && m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr # Avoid fetching terms with discrete variables
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                m.var_type[var_idxs[end]] in [:Bin, :Int] && return false ,expr # Avoid fetching terms with discrete variables
                continue
            end
        end
        # Cofirm detection of patter A and perform store & lifting procedures
        if (length(var_idxs) == 2) && (length(Set(var_idxs)) == 1)
            term_key = [Expr(:ref, :x, var_idxs[1]) for i in 1:2]
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, :MONOMIAL, :*, monomial, basic_monomial_bounds, collect_monomial_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    end

    return false, expr
end

function basic_monomial_bounds(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]
    cnt = 0
    bound = []
    for var in k
        cnt += 1
        var_idx = var.args[2]
        var_bounds = [m.l_var_tight[var_idx], m.u_var_tight[var_idx]]
        if cnt == 1
            bound = copy(var_bounds)
        elseif cnt == 2
            bound = bound * var_bounds'
        else
            bound = diag(bound) * var_bounds'
        end
    end
    if minimum(bound) > m.l_var_tight[lifted_idx] + m.tol
        m.l_var_tight[lifted_idx] = minimum(bound)
    end
    if maximum(bound) < m.u_var_tight[lifted_idx] - m.tol
        m.u_var_tight[lifted_idx] = maximum(bound)
    end

    return
end

function basic_monomial_bounds(m::PODNonlinearModel, nlk::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[nlk][:lifted_var_ref].args[2]
    cnt = 0
    bound = []
    for var in nlk
        cnt += 1
        var_idx = var.args[2]
        var_bounds = [d[var_idx][1], d[var_idx][end]]
        if cnt == 1
            bound = copy(var_bounds)
        elseif cnt == 2
            bound = bound * var_bounds'
        else
            bound = diag(bound) * var_bounds'
        end
    end

    if minimum(bound) > d[lifted_idx][1] + m.tol
        d[lifted_idx][1] = minimum(bound)
    end

    if maximum(bound) < d[lifted_idx][end] - m.tol
        d[lifted_idx][end] = maximum(bound)
    end

    return d
end

function collect_monomial_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)
    for var in k
        @assert isa(var.args[2], Int)
        if var_bowl == nothing
            var.args[2] in m.candidate_disc_vars || push!(m.candidate_disc_vars, var.args[2])
        else
            var.args[2] in var_bowl || push!(var_bowl, var.args[2])
        end
    end
    return
end

"""
    Recognize sin/cos terms: sin(x) / cos(x), where x "should" be continous variables
    # TODO future call detect_TRIGONOMETRIC_term
"""
function detect_sincos_term(expr, constr_id::Int, m::PODNonlinearModel)

    @assert expr.head == :call
    if expr.args[1] in [:sin, :cos]
        # Pattern: sin(a*x) or cos(a*x)
        operator = expr.args[1]
        scalar = 1.0
        var_idxs = []
        for i in 1:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
                continue
            end
            isa(expr.args[i], Symbol) && continue
            (expr.args[i].head == :ref) && isa(expr.args[i].args[2], Int) && push!(var_idxs, expr.args[i].args[2])
            if (expr.args[i].head == :call)
                down_check, linear_lift_var = detect_linear_term(expr.args[i], constr_id, m)
                !down_check && return false, expr
                push!(var_idxs, linear_lift_var.args[2])
                continue
            end
        end
        if length(var_idxs) == 1
            term_key = Dict(:operator=>operator, :scalar=>scalar, :vars=>var_idxs)
            term_key in keys(m.nonconvex_terms) || store_nonconvex_term(m, term_key, var_idxs, term_key[:operator], term_key[:operator], sincos, basic_sincos_bounds, collect_sincos_discvar)
            return true, lift_nonconvex_term(m, term_key, constr_id, scalar)
        end
    end

    return false, expr
end

function basic_sincos_bounds(m::PODNonlinearModel, k::Any)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]
    m.l_var_tight[lifted_idx] = -1  # TODO can be improved
    m.u_var_tight[lifted_idx] = 1

    return
end

function basic_sincos_bounds(m::PODNonlinearModel, k::Any, d::Dict)

    lifted_idx = m.nonconvex_terms[k][:lifted_var_ref].args[2]
    d[lifted_idx][1] = -1  # TODO can be improved
    d[lifted_idx][end] = 1

    return d
end


function collect_sincos_discvar(m::PODNonlinearModel, k::Any; var_bowl=nothing)
    for var in m.nonconvex_terms[k][:var_idxs]
        @assert isa(var, Int)
        if var_bowl == nothing
            var in m.candidate_disc_vars || push!(m.candidate_disc_vars, var)
        else
            var in var_bowl || push!(var_bowl, var)
        end
    end
    return
end

"""
    Recognize convex constraints
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

        # [BUG FIX | FORK]
        power_check = [((i > 1.0) && mod(i, 2) == 0) for i in power_bin]    # Special case for linear constraints with @NLconstraint
        isempty(power_check) && return false

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

        # If we really cannot figure out what is going on...
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
        m.constr_structure[idx] = :convex

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
        m.obj_structure = :convex

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
