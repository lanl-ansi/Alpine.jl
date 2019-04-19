"""
    populate_redundant_constraint_ids!(model::MOI.AbstractOptimizer)

This function identifies the redundant constraints in the model, if any 
"""
function populate_redundant_constraint_ids!(model::MOI.AbstractOptimizer)

    # populate linear redundant constraints
    offset = linear_offset(model)
    num_constraints = model.inner.num_linear_constraints 
    for i in 1:num_constraints 
        func, set = model.linear_constraints[i]
        if isa(set, MOI.LessThan{Float64})
            ub = set.upper 
            new_ub = model.inner.constraint_bound_info[i+offset].hi
            (new_ub < ub) && (push!(model.inner.redundant_constraint_ids.linear_constraint_ids, i))
        elseif isa(set, MOI.GreaterThan{Float64})
            lb = set.lower 
            new_lb = model.inner.constraint_bound_info[i+offset].lo 
            (new_lb > lb) && (push!(model.inner.redundant_constraint_ids.linear_constraint_ids, i))
        else 
            continue
        end 
    end 

    # populate quadratic redundant constraints 
    offset = quadratic_offset(model)
    num_constraints = model.inner.num_quadratic_constraints 
    for i in 1:num_constraints 
        func, set = model.quadratic_constraints[i]
        if isa(set, MOI.LessThan{Float64})
            ub = set.upper 
            new_ub = model.inner.constraint_bound_info[i+offset].hi
            (new_ub < ub) && (push!(model.inner.redundant_constraint_ids.quadratic_constraint_ids, i))
        elseif isa(set, MOI.GreaterThan{Float64})
            lb = set.lower 
            new_lb = model.inner.constraint_bound_info[i+offset].lo 
            (new_lb > lb) && (push!(model.inner.redundant_constraint_ids.quadratic_constraint_ids, i))
        else 
            continue 
        end 
    end

    # handle SOCs and RSOCs later 
    offset = soc_offset(model)
    num_constraints = model.inner.num_soc_constraints
    for i in 1:num_constraints 
        continue
    end 

    offset = rsoc_offset(model)
    num_constraints = model.inner.num_rsoc_constraints
    for i in 1:num_constraints 
        continue
    end 

    # populate nl redundant constraints 
    offset = nl_offset(model)
    num_constraints = model.inner.num_nl_constraints
    for i in 1:num_constraints 
        lb = model.nlp_data.constraint_bounds[i].lower
        ub = model.nlp_data.constraint_bounds[i].upper
        new_lb = model.inner.constraint_bound_info[i+offset].lo 
        new_ub = model.inner.constraint_bound_info[i+offset].hi 

        (!isinf(lb) && new_lb > lb) && (push!(model.inner.redundant_constraint_ids.nl_constraint_ids, i))
        (!isinf(ub) && new_ub < ub) && (push!(model.inner.redundant_constraint_ids.nl_constraint_ids, i))
    end 

    return

end