"""
Forward propagation for linear constraints 
""" 
function fbbt_forward_linear_constraints!(model::MOI.AbstractOptimizer)

    for i in 1:model.inner.num_linear_le_constraints 
        func, set = model.linear_le_constraints[i]
        offset = linear_le_offset(model)
        computed_bound_info = func.constant..func.constant 
        for term in func.terms 
            var_id = term.variable_index.value 
            coeff = term.coefficient 
            computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
        end
        model.inner.constraint_bound_info[offset + i] = 
            intersect(model.inner.constraint_bound_info[offset + i], computed_bound_info)
        if isempty(model.inner.constraint_bound_info[offset + i])
            model.inner.status = Status() 
            model.inner.status.alpine_status = :infeasible
        end 
    end 

    for i in 1:model.inner.num_linear_ge_constraints 
        func, set = model.linear_ge_constraints[i]
        offset = linear_ge_offset(model)
        computed_bound_info = func.constant..func.constant 
        for term in func.terms 
            var_id = term.variable_index.value 
            coeff = term.coefficient  
            computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
        end 
        model.inner.constraint_bound_info[offset + i] = 
            intersect(model.inner.constraint_bound_info[offset + i], computed_bound_info) 
        if isempty(model.inner.constraint_bound_info[offset + i])
            (isa(model.inner.status, Nothing)) && (model.inner.status = Status())
            model.inner.status.alpine_status = :infeasible
        end
    end

    for i in 1:model.inner.num_linear_eq_constraints
        func, set = model.linear_eq_constraints[i]
        offset = linear_eq_offset(model)
        computed_bound_info = func.constant..func.constant 
        for term in func.terms 
            var_id = term.variable_index.value 
            coeff = term.coefficient  
            computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
        end 
        model.inner.constraint_bound_info[offset + i] = 
            intersect(model.inner.constraint_bound_info[offset + i], computed_bound_info) 
        if isempty(model.inner.constraint_bound_info[offset + i])
            (isa(model.inner.status, Nothing)) && (model.inner.status = Status())
            model.inner.status.alpine_status = :infeasible
        end
    end

    return

end 


"""
Backward propagation for linear constraints 
"""
function fbbt_backward_linear_constraints!(model::MOI.AbstractOptimizer)

    # Reference: https://doi.org/10.1080/10556788.2017.1335312 (section 2.1)
    for i in 1:model.inner.num_linear_le_constraints
        func, set = model.linear_le_constraints[i]
        offset = linear_le_offset(model)
        lb = model.inner.constraint_bound_info[offset + i].lo 
        ub = model.inner.constraint_bound_info[offset + i].hi
        for term in func.terms
            id = term.variable_index.value 
            var_lb = lb 
            var_ub = ub 

            for other_term in func.terms 
                other_term_id = other_term.variable_index.value 
                (other_term_id == id) && (continue)
                other_term_coeff = other_term.coefficient 
                other_term_lb = model.inner.variable_bound_original[other_term_id].lo 
                other_term_ub = model.inner.variable_bound_original[other_term_id].hi 
                var_lb -= max(other_term_coeff * other_term_lb, other_term_coeff * other_term_ub)
                var_ub -= min(other_term_coeff * other_term_lb, other_term_coeff * other_term_ub)
            end 

            var_ub = var_ub/term.coefficient
            var_lb = var_lb/term.coefficient

            if term.coefficient > 0.0 
                model.inner.variable_bound_tightened[id] = 
                    intersect(model.inner.variable_bound_tightened[id], var_lb..var_ub)
            end 
            if term.coefficient < 0.0 
                model.inner.variable_bound_tightened[id] = 
                    intersect(model.inner.variable_bound_tightened[id], var_ub..var_lb)
            end
        end 
    end

    for i in 1:model.inner.num_linear_ge_constraints
        func, set = model.linear_ge_constraints[i]
        offset = linear_ge_offset(model)
        lb = model.inner.constraint_bound_info[offset + i].lo 
        ub = model.inner.constraint_bound_info[offset + i].hi
        for term in func.terms
            id = term.variable_index.value 
            var_lb = lb 
            var_ub = ub 

            for other_term in func.terms 
                other_term_id = other_term.variable_index.value 
                (other_term_id == id) && (continue)
                other_term_coeff = other_term.coefficient 
                other_term_lb = model.inner.variable_bound_original[other_term_id].lo 
                other_term_ub = model.inner.variable_bound_original[other_term_id].hi 
                var_lb -= max(other_term_coeff * other_term_lb, other_term_coeff * other_term_ub)
                var_ub -= min(other_term_coeff * other_term_lb, other_term_coeff * other_term_ub)
            end 

            var_ub = var_ub/term.coefficient
            var_lb = var_lb/term.coefficient

            if term.coefficient > 0.0 
                model.inner.variable_bound_tightened[id] = 
                    intersect(model.inner.variable_bound_tightened[id], var_lb..var_ub)
            end 
            if term.coefficient < 0.0 
                model.inner.variable_bound_tightened[id] = 
                    intersect(model.inner.variable_bound_tightened[id], var_ub..var_lb)
            end
        end 
    end

    for i in 1:model.inner.num_linear_eq_constraints
        func, set = model.linear_eq_constraints[i]
        offset = linear_eq_offset(model)
        lb = model.inner.constraint_bound_info[offset + i].lo 
        ub = model.inner.constraint_bound_info[offset + i].hi
        for term in func.terms
            id = term.variable_index.value 
            var_lb = lb 
            var_ub = ub 

            for other_term in func.terms 
                other_term_id = other_term.variable_index.value 
                (other_term_id == id) && (continue)
                other_term_coeff = other_term.coefficient 
                other_term_lb = model.inner.variable_bound_original[other_term_id].lo 
                other_term_ub = model.inner.variable_bound_original[other_term_id].hi 
                var_lb -= max(other_term_coeff * other_term_lb, other_term_coeff * other_term_ub)
                var_ub -= min(other_term_coeff * other_term_lb, other_term_coeff * other_term_ub)
            end 

            var_ub = var_ub/term.coefficient
            var_lb = var_lb/term.coefficient

            if term.coefficient > 0.0 
                model.inner.variable_bound_tightened[id] = 
                    intersect(model.inner.variable_bound_tightened[id], var_lb..var_ub)
            end 
            if term.coefficient < 0.0 
                model.inner.variable_bound_tightened[id] = 
                    intersect(model.inner.variable_bound_tightened[id], var_ub..var_lb)
            end
        end 
    end

    return
end 

"""
Check if problem is infeasible
"""
function is_infeasible(model::MOI.AbstractOptimizer)::Bool 
    if ~isa(model.inner.status, Nothing) && model.inner.status.alpine_status == :infeasible
        return true
    end

    return false
end 