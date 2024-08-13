"""
process_expr(expr; kwargs...)

High-level wrapper for processing expression with sub-tree operators
"""
function process_expr(m::Optimizer)
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
function expr_initialization(m::Optimizer)

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
function expr_preprocess(m::Optimizer)
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
function expr_parsing(m::Optimizer)

    # Throw an error if the objective expression has fractional exponents
    expr_is_fractional_exponent(m.bounding_obj_expr_mip)
    is_structural = expr_constr_parsing(m.bounding_obj_expr_mip, m)

    if !is_structural
        m.bounding_obj_expr_mip = expr_term_parsing(m.bounding_obj_expr_mip, 0, m)
        m.obj_structure = :generic_linear
    end
    (get_option(m, :log_level) > 199) && println("[OBJ] $(m.obj_expr_orig)")

    for i in 1:m.num_constr_orig
        expr = m.bounding_constr_expr_mip[i]

        # Throw an error if the constraint expression has fractional exponents
        if expr.args[1] in [:(<=), :(>=), :(==)]
            expr_is_fractional_exponent(expr.args[2])
        end

        is_structural = expr_constr_parsing(expr, m, i)
        if !is_structural
            m.bounding_constr_expr_mip[i] = expr_term_parsing(expr, i, m)
            m.constr_structure[i] = :generic_linear
        end
        (get_option(m, :log_level) > 199) && println("[CONSTR] $(m.constr_expr_orig[i])")
    end

    return
end

"""
STEP 4: convert the parsed expressions into affine-based function that can be used for adding JuMP constraints
"""
function expr_conversion(m::Optimizer)
    if m.obj_structure == :generic_linear
        m.bounding_obj_mip = expr_linear_to_affine(m.bounding_obj_expr_mip)
        m.obj_structure = :affine
    end
    get_option(m, :log_level) > 199 && println("type :: ", m.obj_structure)
    get_option(m, :log_level) > 199 && println("lifted ::", m.bounding_obj_expr_mip)
    get_option(m, :log_level) > 199 && println("coeffs ::", m.bounding_obj_mip[:coefs])
    get_option(m, :log_level) > 199 && println("vars ::", m.bounding_obj_mip[:vars])
    get_option(m, :log_level) > 199 && println("sense ::", m.bounding_obj_mip[:sense])
    get_option(m, :log_level) > 199 && println("rhs ::", m.bounding_obj_mip[:rhs])
    get_option(m, :log_level) > 199 && println("----------------")

    for i in 1:m.num_constr_orig
        if m.constr_structure[i] == :generic_linear
            m.bounding_constr_mip[i] =
                expr_linear_to_affine(m.bounding_constr_expr_mip[i])
            m.constr_structure[i] = :affine
        end
        get_option(m, :log_level) > 199 && println("type :: ", m.constr_structure[i])
        get_option(m, :log_level) > 199 &&
            println("lifted ::", m.bounding_constr_expr_mip[i])
        get_option(m, :log_level) > 199 &&
            println("coeffs ::", m.bounding_constr_mip[i][:coefs])
        get_option(m, :log_level) > 199 &&
            println("vars ::", m.bounding_constr_mip[i][:vars])
        get_option(m, :log_level) > 199 &&
            println("sense ::", m.bounding_constr_mip[i][:sense])
        get_option(m, :log_level) > 199 &&
            println("rhs ::", m.bounding_constr_mip[i][:rhs])
        get_option(m, :log_level) > 199 && println("----------------")
    end

    return
end

"""
STEP 5: collect information as needed for handy operations in the algorithm section
"""
function expr_finalized(m::Optimizer)
    collect_nonconvex_vars(m)
    m.candidate_disc_vars = sort(m.candidate_disc_vars)
    m.num_var_linear_mip = length(m.linear_terms)
    m.num_var_nonlinear_mip = length(m.nonconvex_terms)
    m.num_constr_convex = length([i for i in m.constr_structure if i == :convex])

    return m
end

function collect_nonconvex_vars(m::Optimizer)

    # Walk through all nonconvex terms
    for i in keys(m.nonconvex_terms)
        m.nonconvex_terms[i][:discvar_collector](m, i)
    end

    # TODO : reconsider how to structure this
    # NOTE : Walk through all integer variables that didn't appear in any nonconvex terms
    for i in 1:m.num_var_orig
        if !(i in m.candidate_disc_vars) && m.var_type[i] == :Int
            push!(m.candidate_disc_vars, i)
        end
    end

    return
end

function isa_variable_index(expr::Expr)
    return length(expr.args == 2) && expr.args[2] == :(MathOptInterface.VariableIndex)
end
get_index(expr::Expr) = expr.args[2]

function expr_strip_const(expr, subs = [], rhs = 0.0)
    exhaust_const = [
        !(expr.args[1] in [:+, :-]) ||
        !(isa(expr.args[i], Float64) || isa(expr.args[i], Int)) for
        i in 2:length(expr.args)
    ]
    if prod(exhaust_const)
        push!(subs, expr)
        return subs, rhs
    end

    for i in 2:length(expr.args)
        if (
            isa(expr.args[i], Float64) ||
            isa(expr.args[i], Int) ||
            isa(expr.args[i], Symbol)
        )
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

"""
This utility function builds a constraint reference by repeating one
operator with a vector variable references.
For example,
input => y, x[1,2,3], :+
output => (y = x[1] + x[2] + x[3])::Expr
"""
function build_constr_block(y_idx::Int, var_idxs::Vector, operator::Symbol)
    expr_l_block = Expr(:call, :(==), Meta.parse("x[$(y_idx)]"))
    expr_r_block = Expr(:call, operator)
    for j in 1:length(var_idxs)
        push!(expr_r_block.args, Meta.parse("x[$(var_idxs[j])]"))
    end
    push!(expr_l_block.args, expr_r_block)

    return expr_l_block
end

"""
expr_constr_parsing(expr, m::Optimizer)

Recognize structural constraints.
"""
function expr_constr_parsing(expr::Expr, m::Optimizer, idx::Int = 0)

    # First process user-defined structures in-cases of over-ride
    # for i in 1:length(get_option(m, :constr_patterns))
    #    is_structural = get_option(m, :constr_patterns)[i](expr, m, idx)
    #    return
    # end

    isa(expr, Number) && return false
    # Recognize built-in special structural pattern
    if get_option(m, :recognize_convex)
        is_convex = resolve_convex_constr(expr, m, idx)
        is_convex && return true
    end

    # More structural patterns go here
    return false
end

function expr_is_axn(expr, scalar = 1.0, var_idxs = [], power = []; N = nothing)
    expr.args[1] in [:*, :^] || return nothing, nothing, nothing         # Limited Area

    if expr.args[1] == :*
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
            elseif (expr.args[i].head == :ref)
                @assert isa(expr.args[i].args[2], Int)
                push!(var_idxs, expr.args[i])
                push!(power, 1)
            elseif (expr.args[i].head == :call)
                scalar, var_idxs, power =
                    expr_is_axn(expr.args[i], scalar, var_idxs, power)
                scalar === nothing && return nothing, nothing, nothing
            end
        end
    elseif expr.args[1] == :^
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                push!(power, expr.args[i])
                continue
            elseif (expr.args[i].head == :ref)
                push!(var_idxs, expr.args[i])
                continue
            elseif (expr.args[i].head == :call)
                return nothing, nothing, nothing
            end
        end
    end

    # If the user wants a specific N
    !(N === nothing) && !(length(var_idxs) == N) && return nothing, nothing, nothing
    (var_idxs === nothing) &&
        (scalar === nothing) &&
        (power === nothing) &&
        return nothing, nothing, nothing # Unrecognized sub-structure

    @assert length(var_idxs) == length(power)
    return scalar, var_idxs, power
end

"""
This function takes a constraint/objective expression and converts it into a affine expression data structure
Use the function to traverse linear expressions traverse_expr_linear_to_affine()
"""
function expr_linear_to_affine(expr)

    # The input should follow :(<=, LHS, RHS)
    affdict = Dict()
    isa(expr, Number) && return affdict

    if expr.args[1] in [:(==), :(>=), :(<=)] # For a constraint expression
        @assert isa(expr.args[3], Float64) || isa(expr.args[3], Int)
        @assert isa(expr.args[2], Expr)
        # non are buffer spaces, not used anywhere
        lhscoeff, lhsvars, rhs, non, non = traverse_expr_linear_to_affine(expr.args[2])
        rhs = -rhs + expr.args[3]
        affdict[:sense] = expr.args[1]

    elseif expr.head == :ref  # For single variable objective expression
        lhscoeff = [1.0]
        lhsvars = [expr]
        rhs = 0
        affdict[:sense] = nothing

    else # For an objective expression
        lhscoeff, lhsvars, rhs, non, non = traverse_expr_linear_to_affine(expr)
        affdict[:sense] = nothing
    end

    affdict[:coefs] = lhscoeff
    affdict[:vars] = lhsvars
    affdict[:rhs] = rhs
    @assert length(affdict[:coefs]) == length(affdict[:vars])
    affdict[:cnt] = length(affdict[:coefs])

    return affdict
end

"""

traverse_expr_linear_to_affine(expr, lhscoeffs=[], lhsvars=[], rhs=0.0, bufferVal=0.0, bufferVar=nothing, sign=1.0, level=0)

This function traverses the left hand side tree to collect affine terms.
Updated status : possible to handle (x-(x+y(t-z))) cases where signs are handled properly
"""
function traverse_expr_linear_to_affine(
    expr,
    lhscoeffs = [],
    lhsvars = [],
    rhs = 0.0,
    bufferVal = nothing,
    bufferVar = nothing,
    sign = 1.0,
    coef = 1.0,
    level = 0,
)

    # reversor = Dict(true => -1.0, false => 1.0)
    function sign_convertor(subexpr, pos)
        if length(subexpr.args) == 2 && subexpr.args[1] == :-
            return -1.0
        elseif length(subexpr.args) > 2 && subexpr.args[1] == :- && pos > 2
            return -1.0
        end
        return 1.0
    end

    if isa(expr, Number) # Capture any coefficients or right hand side
        (bufferVal !== nothing) ? bufferVal *= expr : bufferVal = expr * coef
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:+, :-]    # TODO: what is this condition?
        if bufferVal !== nothing && bufferVar !== nothing
            push!(lhscoeffs, bufferVal)
            push!(lhsvars, bufferVar)
            bufferVal = 0.0
            bufferVar = nothing
        end
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:*]
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:(<=), :(==), :(>=)]
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr.head == :ref
        bufferVar = expr
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    end

    # PATCH : Special Structure Recognition
    start_pos = 1
    if (expr.args[1] == :*) && (length(expr.args) == 3)
        if (isa(expr.args[2], Float64) || isa(expr.args[2], Int)) &&
           (expr.args[3].head == :call)
            (coef != 0.0) ? coef = expr.args[2] * coef : coef = expr.args[2]# Patch
            # coef = expr.args[2]
            start_pos = 3
            # @warn "Special expression structure detected [*, coef, :call, ...]. Currently using a beta fix..."
        end
    end

    for i in start_pos:length(expr.args)
        lhscoeff, lhsvars, rhs, bufferVal, bufferVar = traverse_expr_linear_to_affine(
            expr.args[i],
            lhscoeffs,
            lhsvars,
            rhs,
            bufferVal,
            bufferVar,
            sign * sign_convertor(expr, i),
            coef,
            level + 1,
        )
        if expr.args[1] in [:+, :-]  # Term segmentation [:-, :+], see this and wrap-up the current (linear) term
            if bufferVal !== nothing && bufferVar !== nothing  # (sign) * (coef) * (var) => linear term
                push!(lhscoeffs, sign * sign_convertor(expr, i) * bufferVal)
                push!(lhsvars, bufferVar)
                bufferVal = nothing
                bufferVar = nothing
            end
            if bufferVal !== nothing && bufferVar === nothing  # (sign) * (coef) => right-hand-side term
                rhs += sign * sign_convertor(expr, i) * bufferVal
                bufferVal = nothing
            end
            if bufferVal === nothing && bufferVar !== nothing && expr.args[1] == :+
                push!(lhscoeffs, sign * 1.0 * coef)
                push!(lhsvars, bufferVar)
                bufferVar = nothing
            end
            if bufferVal === nothing && bufferVar !== nothing && expr.args[1] == :-
                push!(lhscoeffs, sign * sign_convertor(expr, i) * coef)
                push!(lhsvars, bufferVar)
                bufferVar = nothing
            end
        elseif expr.args[1] in [:(<=), :(==), :(>=)]
            rhs = expr.args[end]
        end
    end

    if level == 0
        if bufferVal !== nothing && bufferVar !== nothing
            push!(lhscoeffs, bufferVal)
            push!(lhsvars, bufferVar)
            bufferVal = nothing
            bufferVar = nothing
        end
    end

    return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
end

"""
This function is treats the sign in expression trees by cleaning the following structure:
args
::+ || ::-
::Call || ::ref
By separating the structure with some dummy treatments
"""
function expr_resolve_sign(expr, level = 0; kwargs...)
    isa(expr, Number) && return
    resolver = Dict(:- => -1, :+ => 1)
    for i in 2:length(expr.args)
        if !isa(expr.args[i], Float64) && !isa(expr.args[i], Int) # Skip the coefficients
            if length(expr.args[i].args) == 2 && expr.args[i].head == :call
                if expr.args[i].args[1] == :- # Treatment for :(-x), replace with :(x*-1)
                    if expr.args[1] in [:*] # Only perform the treatment when connected with *
                        push!(expr.args, resolver[expr.args[i].args[1]])
                        expr.args[i] = expr.args[i].args[2]
                    elseif expr.args[1] in [:+, :-, :(==), :(<=), :(>=)]
                        expr.args[i] = expr.args[i]
                    else
                        error("Unexpected operator $(expr.args[i]) during resolving sign")
                    end
                elseif expr.args[i].args[1] == :+ # Treatment for :(+x) replace with :x
                    expr.args[i] = expr.args[i].args[2]
                end
            elseif expr.args[i].head == :call
                expr_resolve_sign(expr.args[i], level + 1)
            end
        end
    end

    return
end

"""
Recursively pre-process the expression by treating the coefficients
TODO: this function requires a lot of refining.
Most issues can be caused by this function.
"""
function expr_flatten(expr, level = 0; kwargs...)
    isa(expr, Number) && return
    if level > 0  # No trivial constraint is allowed "3>5"
        flat = expr_arrangeargs(expr.args)
        if isa(flat, Number)
            return flat
        else
            expr.args = flat
        end
    end

    for i in 2:length(expr.args)
        if isa(expr.args[i], Expr) && expr.args[i].head == :call
            expr.args[i] = expr_flatten(expr.args[i], level + 1)
        end
    end

    if level > 0  #Root level process, no additional processes
        expr.args = expr_arrangeargs(expr.args)
    end

    return expr
end

"""
Re-arrange children by their type.
Only consider 3 types: coefficients(Int64, Float64), :ref, :calls
Sequence arrangement is dependent on operator
"""
function expr_arrangeargs(args::Array; kwargs...)
    # A mapping used for operator treatments
    reverter = Dict(true => :-, false => :+)
    operator_coefficient_base =
        Dict{Symbol,Float64}(:+ => 0.0, :- => 0.0, :/ => 1.0, :* => 1.0, :^ => 1.0)

    # Treatment => Flatten pure number sub-tree
    allFloats = [i for i in 2:length(args) if isa(args[i], Float64)]
    if length(allFloats) == length(args) - 1
        val = args[2]
        for i in 3:length(args)
            val = eval(args[1])(val, args[i])
        end
        return val
    end

    if args[1] in [:^]
        return args
    elseif args[1] in [:/]
        # error("Alpine does not currently support `$(args[1])` operator")
        if (typeof(args[3]) == Float64) || (typeof(args[3]) == Int64)
            args = expr_resolve_const_denominator(args)
        else
            error(
                "Alpine does not currently support `$(args[1])` operator with a variable in the denominator",
            )
        end
    elseif args[1] in [:*, :+, :-]
        # Such generic operators can be handled
    else
        error(
            "Alpine does not currently support `$(args[1])` operator acting on a variable",
        )
    end

    if args[1] in [:*, :+, :-]
        val = operator_coefficient_base[args[1]] # initialize
        exprs = []
        refs = []

        valinit = false  # what is first children : needs consideration with :(-)
        refinit = false
        callinit = false

        for i in 2:length(args)
            if isa(args[i], Float64) || isa(args[i], Int)
                valinit += (i == 2)
                val = eval(args[1])(val, eval(reverter[(i==2)&&(args[1]==:-)])(args[i])) # Revert the first sigh if (-)
            elseif typeof(args[i]) == Expr
                if args[i].head == :ref
                    refinit += (i == 2)
                    push!(refs, args[i])
                elseif args[i].head == :call
                    callinit += (i == 2)
                    push!(exprs, args[i])
                else
                    error("Unkown expression head")
                end
            else
                error("Unkown argument type")
            end
        end

        # Treatment for 0.0 coefficients
        if val == operator_coefficient_base[args[1]] && val == 0.0
            val = []
        end

        # Re-arrange children :: without actually doing anything on them
        if args[1] in [:*]
            return [args[1]; val; refs; exprs]
        elseif args[1] in [:+, :-]
            if Bool(valinit)
                return [args[1]; val; refs; exprs]
            elseif Bool(refinit)
                return [args[1]; refs; val; exprs]
            elseif Bool(callinit)
                return [args[1]; exprs; refs; val]
            else
                error("Unexpected condition encountered...")
            end
        end
    else
        error("Unsupported expression arguments $args")
    end

    return
end

"""
Check if it is a sample of `+()`
"""
expr_is_emptysum(expr) = expr == :(+())

"""
Check if a sub-tree is a constant or not
"""
function expr_resolve_const(expr)
    isa(expr, Number) && return
    for i in 1:length(expr.args)
        if isa(expr.args[i], Number) || isa(expr.args[i], Symbol)
            continue
        elseif expr.args[i].head == :call
            if expr_is_emptysum(expr.args[i])
                expr.args[i] = :(0)
            end
            (expr_isconst(expr.args[i])) && (expr.args[i] = eval(expr.args[i]))
            if isa(expr.args[i], Number) || isa(expr.args[i], Symbol)
                continue
            else
                expr_resolve_const(expr.args[i])
            end
        end
    end

    return
end

"""
If the expression's fraction has a constant denominator, then replace the expression as product of a constant and the expression
"""
function expr_resolve_const_denominator(args)
    @assert args[1] == :/
    @assert length(args) == 3
    resolvable = true

    for i in 3:length(args)
        resolvable *= expr_isconst(args[i])
    end

    if resolvable
        for i in 3:length(args)
            args[i] = 1 / eval(args[i])
        end
        args[1] = :*
    end

    return args
end

"""
Check if a sub-tree(:call) is totally composed of constant values
"""
function expr_isconst(expr)
    (isa(expr, Number) || isa(expr, Symbol)) && return true

    (expr.head == :ref) && return false
    const_tree = true
    for i in 1:length(expr.args)
        if isa(expr.args[i], Number) || isa(expr.args[i], Symbol)
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
Check if a sub-tree(:call) is contains any non-integer exponent values
"""
function expr_is_fractional_exponent(expr)
    if expr.head == :call
        if length(expr.args) == 3 && expr.args[1] == :^
            if !(isinteger(expr.args[3])) || !(expr.args[3] >= 0)
                error(
                    "Alpine currently supports `^` operator in constraints and/or objective with only positive integer exponents",
                )
            end
        end
        for i in 1:length(expr.args)
            if typeof(expr.args[i]) == Expr
                expr_is_fractional_exponent(expr.args[i]) # Recursively search for fractional exponents
            end
        end
    end
end

# Returns true if the expression is a constant, linear or affine
function expr_isaffine(expr)
    expr_isconst(expr) && return true
    expr.head == :ref && return true

    is_affine = false
    if expr.head == :call
        k = 0
        for i in 1:length(expr.args)
            if isa(expr.args[i], Number)
                k += 1
                continue
            end
            if isa(expr.args[i], Symbol)
                (expr.args[i] in [:+, :-]) && (k += 1)
                continue
            end
            if expr.args[i].head == :ref
                k += 1
            elseif expr.args[i].head == :call
                status = expr_isaffine(expr.args[i])
                status && (k += 1)
            end
        end
        (k == length(expr.args)) && (is_affine = true)
    end
    return is_affine
end

"""
Converts ((a_1*x[1])^2 + (a_2*x[2])^2 + ... + (a_n*x[n])^2) to (a_1^2*x[1]^2 + a_2^2*x[2]^2 + ... + a_n^2*x[n]^2)
Signs in the summation can be +/-
Note: This function does not support terms of type (a*(x[1] + x[2]))^2 yet.
"""
function expr_isolate_const(expr)
    expr_isaffine(expr) && return expr

    # Check nonlinear expressions
    if (expr.head == :call && expr.args[1] in [:+, :-])
        expr_array = Any[]
        for i in 2:length(expr.args)
            ind = 0

            # Handle negative sign in the first term
            if (!isa(expr.args[i], Number)) &&
               (expr.args[i].args[1] == :-) &&
               (length(expr.args[i].args) == 2)
                expr_i = expr.args[i].args[2]
                ind = 1
            else
                expr_i = expr.args[i]
            end

            # Handle constant and linear terms
            if (isa(expr_i, Number)) || (expr_i.head == :ref)
                if ((expr.args[1] == :-) && (i > 2)) || (ind == 1)
                    push!(expr_array, :(-$(expr_i)))
                else
                    push!(expr_array, :($(expr_i)))
                end

                # Handle nonlinear terms with coefficients within the exponent
            elseif (
                expr.args[i].head == :call &&
                expr_i.args[1] == :^ &&
                expr_i.args[2].head == :call &&
                isa(expr_i.args[2].args[2], Number) &&
                isa(expr_i.args[3], Number)
            )
                expr_tmp = Expr(:call, :^, expr_i.args[2].args[3], expr_i.args[3])
                if ((expr.args[1] == :-) && (i > 2)) || (ind == 1)
                    push!(
                        expr_array,
                        Expr(:call, :*, -expr_i.args[2].args[2]^expr_i.args[3], expr_tmp),
                    )
                else
                    push!(
                        expr_array,
                        Expr(:call, :*, expr_i.args[2].args[2]^expr_i.args[3], expr_tmp),
                    )
                end

                # Handle no-coefficients case
            elseif expr_i.args[1] == :^ && expr_i.args[2].head == :ref
                if (expr.args[1] == :- && (i > 2)) || (ind == 1)
                    push!(expr_array, Expr(:call, :*, -1, expr_i))
                else
                    push!(expr_array, expr_i)
                end

                # Handle coefficients which are not part of the exponent
            elseif expr_i.args[1] == :* && isa(expr_i.args[2], Number)
                if ((expr.args[1] == :-) && (i > 2)) || (ind == 1)
                    #push!(expr_array, Expr(:call, :*, -expr_i.args[2], expr_i.args[3]))
                    push!(expr_array, Expr(:call, :*, -1, expr_i))
                else
                    push!(expr_array, expr_i)
                end

                # Handle negative sign in the remaining terms
            elseif (expr_i.args[1] in [:+, :-]) && (length(expr.args[i].args) > 2)
                expr_rec = expr_isolate_const(expr_i) #recursion
                push!(expr_array, expr_rec)

                # For any other terms
            else
                if (expr.args[1] == :- && (i > 2)) || (ind == 1)
                    push!(expr_array, Expr(:call, :*, -1, expr_i))
                else
                    push!(expr_array, expr_i)
                end
            end
        end

        # Construct the expression from the array
        if length(expr_array) == 1
            return expr_array[1]
        else
            expr_n = Expr(:call, :+, expr_array[1], expr_array[2])
            if length(expr_array) >= 3
                for i in 3:length(expr_array)
                    expr_n = Expr(:call, :+, expr_n, expr_array[i])
                end
            end
            return (expr_n)
        end

    elseif (
        expr.head == :call &&
        expr.args[1] == :^ &&
        expr.args[2].head == :call &&
        isa(expr.args[2].args[2], Number) &&
        isa(expr.args[3], Number)
    )
        expr_tmp = Expr(:call, :^, expr.args[2].args[3], expr.args[3])
        return (Expr(:call, :*, expr.args[2].args[2]^expr.args[3], expr_tmp))
    else
        return expr
    end
end
