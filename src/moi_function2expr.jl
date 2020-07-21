_variable_index_to_expr(v::MOI.VariableIndex) = Expr(:ref, :x, v)

function _moi_function_to_expr(t::MOI.ScalarAffineTerm)
    return Expr(:call, :*, t.coefficient, _variable_index_to_expr(t.variable_index))
end

function _moi_function_to_expr(f::MOI.ScalarAffineFunction)
    if length(terms(f)) == 1
        return Expr(:call, :+, _moi_function_to_expr(first(terms(f))), constant(f))
    else
        return Expr(:call, 
                    :+, 
                    _moi_function_to_expr(first(terms(f))),
                    _moi_function_to_expr(
                                 MOI.ScalarAffineFunction(terms(f)[2:end], constant(f)))
                   )
    end
end

function _moi_function_to_expr(t::MOI.ScalarQuadraticTerm)
    return Expr(:call, :*, 
                MOI.coefficient(t), 
                Expr(:call, :*, 
                     _variable_index_to_expr(t.variable_index_1),
                     _variable_index_to_expr(t.variable_index_2)
                    )
               )
end

function _moi_function_to_expr(f::MOI.ScalarQuadraticFunction)
    if length(f.quadratic_terms) == 1
        return Expr(:call, :+, 
                    _moi_function_to_expr(first(f.quadratic_terms)),
                    _moi_function_to_expr(MOI.ScalarAffineFunction(f.affine_terms, constant(f)))
                   )
                else
        return Expr(:call, 
                    :+, 
                    _moi_function_to_expr(first(f.quadratic_terms)),
                    _moi_function_to_expr(
                                 MOI.ScalarQuadraticFunction(
                                             f.affine_terms, 
                                             f.quadratic_terms[2:end], 
                                             constant(f)
                                            )
                                 )
                   )
    end
end
