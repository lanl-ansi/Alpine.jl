_variable_index_to_expr(v::MOI.VariableIndex) = Expr(:ref, :x, v.value)

function _add_constant(expr::Expr, constant)
    if !iszero(constant)
        push!(expr.args, constant)
    end
end

function _term_to_expr(t::MOI.ScalarAffineTerm)
    return Expr(:call, :*, t.coefficient, _variable_index_to_expr(t.variable_index))
end

function _add_terms(expr::Expr, terms::Vector{<:MOI.ScalarAffineTerm})
    for term in terms
        push!(expr.args, _term_to_expr(term))
    end
end

function _moi_function_to_expr(f::MOI.ScalarAffineFunction)
    iszero(f) && return 0 # Similar to `JuMP.affToExpr` in JuMP v0.18 and earlier
    expr = Expr(:call, :+)
    _add_terms(expr, f.terms)
    _add_constant(expr, f.constant)
    return expr
end

function _moi_function_to_expr(t::MOI.ScalarQuadraticTerm)
    return Expr(
        :call, :*,
        MOI.coefficient(t),
        _variable_index_to_expr(t.variable_index_1),
        _variable_index_to_expr(t.variable_index_2)
    )
end

function _moi_function_to_expr(f::MOI.ScalarQuadraticFunction)
    iszero(f) && return 0 # Similar to `JuMP.quadToExpr` in JuMP v0.18 and earlier
    expr = Expr(:call, :+)
    _add_terms(expr, f.quadratic_terms)
    _add_terms(expr, f.affine_terms)
    _add_constant(expr, f.constant)
    return expr
end
