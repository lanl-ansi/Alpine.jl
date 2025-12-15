function _constraint_expr(expr, set::MOI.Interval)
    return Expr(:comparison, set.lower, :(<=), expr, :(<=), set.upper)
end

_sense(::MOI.EqualTo) = :(==)
_sense(::MOI.LessThan) = :(<=)
_sense(::MOI.GreaterThan) = :(>=)
function _constraint_expr(expr, set)
    return Expr(:call, _sense(set), expr, MOI.constant(set))
end

_moi_function_to_expr(f::Number) = f

_moi_function_to_expr(v::MOI.VariableIndex) = Expr(:ref, :x, v.value)

function _moi_function_to_expr(t::MOI.ScalarAffineTerm)
    return Expr(:call, :*, t.coefficient, _moi_function_to_expr(t.variable))
end

function _moi_function_to_expr(t::MOI.ScalarQuadraticTerm)
    coef = t.coefficient
    if t.variable_1 == t.variable_2
        coef /= 2
    end
    return Expr(
        :call,
        :*,
        coef,
        _moi_function_to_expr(t.variable_1),
        _moi_function_to_expr(t.variable_2),
    )
end

function _add_constant(expr::Expr, constant)
    if !iszero(constant)
        push!(expr.args, constant)
    end
    return
end

function _add_terms(expr::Expr, terms::Vector)
    for term in terms
        push!(expr.args, _moi_function_to_expr(term))
    end
    return
end

function _moi_function_to_expr(f::MOI.ScalarAffineFunction)
    if iszero(f)
        return 0 # Similar to `JuMP.affToExpr` in JuMP v0.18 and earlier
    end
    expr = Expr(:call, :+)
    _add_terms(expr, f.terms)
    _add_constant(expr, f.constant)
    return expr
end

function _moi_function_to_expr(f::MOI.ScalarQuadraticFunction)
    if iszero(f)
        return 0 # Similar to `JuMP.quadToExpr` in JuMP v0.18 and earlier
    end
    expr = Expr(:call, :+)
    _add_terms(expr, f.quadratic_terms)
    _add_terms(expr, f.affine_terms)
    _add_constant(expr, f.constant)
    return expr
end

function _moi_function_to_expr(f::MOI.ScalarNonlinearFunction)
    return Expr(
        :call,
        f.head,
        _moi_function_to_expr.(f.args)...,
    )
end
