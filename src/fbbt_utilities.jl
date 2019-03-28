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
Forward propogation for quadratic constraints
"""
function fbbt_forward_quadratic_constraints!(model::MOI.AbstractOptimizer)

    for i in 1:model.inner.num_quadratic_le_constraints 
        func, set = model.quadratic_le_constraints[i]
        offset = quadratic_le_offset(model)
        computed_bound_info = func.constant..func.constant 
        for term in func.affine_terms 
            var_id = term.variable_index.value 
            coeff = term.coefficient 
            computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
        end
        for term in func.quadratic_terms 
            var_id_1 = term.variable_index_1.value 
            var_id_2 = term.variable_index_2.value 
            coeff = term.coefficient/2.0
            if var_id_1 == var_id_2 
                computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id_1]^2 
            else 
                computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id_1] * 
                    model.inner.variable_bound_tightened[var_id_2]
            end 
        end 
        model.inner.constraint_bound_info[offset + i] = 
            intersect(model.inner.constraint_bound_info[offset + i], computed_bound_info)
        if isempty(model.inner.constraint_bound_info[offset + i])
            model.inner.status = Status() 
            model.inner.status.alpine_status = :infeasible
        end 
    end 

    for i in 1:model.inner.num_quadratic_ge_constraints 
        func, set = model.quadratic_ge_constraints[i]
        offset = quadratic_ge_offset(model)
        computed_bound_info = func.constant..func.constant 
        for term in func.affine_terms 
            var_id = term.variable_index.value 
            coeff = term.coefficient 
            computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
        end
        for term in func.quadratic_terms 
            var_id_1 = term.variable_index_1.value 
            var_id_2 = term.variable_index_2.value 
            coeff = term.coefficient/2.0
            if var_id_1 == var_id_2 
                computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id_1]^2 
            else 
                computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id_1] * 
                    model.inner.variable_bound_tightened[var_id_2]
            end 
        end 
        model.inner.constraint_bound_info[offset + i] = 
            intersect(model.inner.constraint_bound_info[offset + i], computed_bound_info) 
        if isempty(model.inner.constraint_bound_info[offset + i])
            (isa(model.inner.status, Nothing)) && (model.inner.status = Status())
            model.inner.status.alpine_status = :infeasible
        end
    end

    for i in 1:model.inner.num_quadratic_eq_constraints
        func, set = model.quadratic_eq_constraints[i]
        offset = quadratic_eq_offset(model)
        computed_bound_info = func.constant..func.constant 
        for term in func.affine_terms 
            var_id = term.variable_index.value 
            coeff = term.coefficient 
            computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
        end
        for term in func.quadratic_terms 
            var_id_1 = term.variable_index_1.value 
            var_id_2 = term.variable_index_2.value 
            coeff = term.coefficient/2.0
            if var_id_1 == var_id_2 
                computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id_1]^2 
            else 
                computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id_1] * 
                    model.inner.variable_bound_tightened[var_id_2]
            end 
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
Forward propagation for DAG; Reference: https://doi.org/10.1007/s11590-019-01396-y  
""" 
function fbbt_forward_dag!(model::MOI.AbstractOptimizer)

    dag = model.inner.expression_graph 
    dag_lookup = model.inner.dag_lookup 
    variable_bound = model.inner.variable_bound_tightened

    max_depth = dag.max_depth 

    for dag_vertex in dag.vertices[0]
        if dag_vertex.vertex_type == :variable 
            id = dag_vertex.vertex.args[2].value
            dag_vertex.interval = intersect(dag_vertex.interval, variable_bound[id])
        end 
    end 

    for d in 1:max_depth 
        for dag_vertex in dag.vertices[d]
            operation = dag_vertex.vertex.args[1]
            children = dag_vertex.children 
            children_interval = Any[]
            for child in children
                depth, position = dag_lookup[child]
                push!(children_interval, dag.vertices[depth][position].interval)
            end 
            if length(children_interval) == 1
                new_interval = eval(Expr(:call, operation, children_interval[1]))
                dag_vertex.interval = intersect(dag_vertex.interval, new_interval)
            elseif length(children_interval) == 2 
                new_interval = eval(Expr(:call, operation, children_interval[1], children_interval[2]))
                dag_vertex.interval = intersect(dag_vertex.interval, new_interval)
            else 
                new_expr = Expr(:call, operation, children_interval[1], children_interval[2])
                for i in 3:length(children_interval)
                    push!(new_expr.args, children_interval[i])
                end 
                new_interval = eval(new_expr)
                dag_vertex.interval = intersect(dag_vertex.interval, new_interval)
            end 
        end 
    end

    return
end 

"""
Forward propogation for nonlinear constraints 
"""
function fbbt_forward_nl_constraints!(model::MOI.AbstractOptimizer)
    fnames = fieldnames(NLFunction)
    dag_lookup = model.inner.dag_lookup
    dag = model.inner.expression_graph
    for i in 1:model.inner.num_nlp_constraints 
        nl_function = model.inner.nl_function[i]
        offset = nlp_offset(model)
        computed_bound_info = 0..0
        for fname in fnames
            f_value = getfield(nl_function, fname)
            (isa(f_value, Nothing)) && (continue)
            for alpine_expression in f_value
                coeff, expr = alpine_expression.expression
                if fname == :constant_part 
                    computed_bound_info += coeff * expr..expr  
                elseif fname == :linear_part 
                    var_id = expr.args[2].value
                    computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
                else 
                    depth, position = depth, position = dag_lookup[expr]
                    computed_bound_info += coeff * dag.vertices[depth][position].interval
                end 
            end 
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
Backward propogation for DAG 
"""
function fbbt_backward_dag!(model::MOI.AbstractOptimizer)
    
    dag = model.inner.expression_graph 
    dag_lookup = model.inner.dag_lookup 
    
    max_depth = dag.max_depth
    
    for d in max_depth:-1:2
        for dag_vertex in dag.vertices[d]
            operation = dag_vertex.vertex.args[1]
            current_interval = dag_vertex.interval
            children = dag_vertex.children 
            children_id = Tuple{Int, Int}[]
            children_interval = Interval{Float64}[] 
            for child in children 
                depth, position = dag_lookup[child]
                children_interval[(depth, position)] = dag.vertices[depth][position].interval
            end 
            # propagate_intervals!(operation, current_interval, children_interval)
        end 

    end 

    return
    
end 

"""
Propagate intervals backward given an operation, current_interval and intervals for each child 
Reference: https://edoc.hu-berlin.de/bitstream/handle/18452/17356/vigerske.pdf?sequence=1&isAllowed=y (Page 175)

function propagate_intervals!(operation::Symbol, current_interval::Interval{Float64}, children_interval::Vector{Interval{Float64}})

    if operation == :- && length(children_interval) == 1
        children_interval[1] = -current_interval
    elseif operation == :- && length(children_interval) == 2 
        current_interval, children_interval[1], children_interval[2] = minus_rev(current_interval, children_interval[1], children_interval[2])
    elseif operation == :+ && length(childrent_interval) == 2
        current_interval, children_interval[1], children_interval[2] = plus_rev(current_interval, children_interval[1], children_interval[2])
    elseif operation == :* && length(childrent_interval) == 2
        current_interval, children_interval[1], children_interval[2] = mul_rev(current_interval, children_interval[1], children_interval[2])
    elseif operation == :/ 
        current_interval, children_interval[1], children_interval[2] = div_rev(current_interval, children_interval[1], children_interval[2])
    elseif operation == :exp 
        lo = max(current_interval.lo, 0.0)
        hi = current_interval.hi
        children_interval[1] = log(lo..hi) 
    elseif operation == :log 
        children_interval[1] = exp(current_interval) 
    elseif operation == :sqrt 
        current_interval, children_interval[1] = sqrt_rev(current_interval, children_interval[1])
    elseif operation == :^
        is_integer = try isa(Int(children[2]), Int) catch y false end 
        (is_integer && Int(children[2]) == 2) && (current_interval, children_interval[1] = sqr_rev(current_interval, children_interval[1]))
        (is_integer && Int(children[2]) != 2) && (current_interval, children_interval[1], n = power_rev(current_interval, children_interval[1], Int(children[2])))
        if !is_integer
            lo = max(current_interval.lo, 0) 
            hi = current_interval.hi 
            a = children_interval[2].lo
            children_interval[1] = lo..hi^(1/a)
        end
    end 
    return 
end 
"""

"""
Check if problem is infeasible
"""
function is_infeasible(model::MOI.AbstractOptimizer)::Bool 
    if ~isa(model.inner.status, Nothing) && model.inner.status.alpine_status == :infeasible
        return true
    end

    return false
end 

"""
for d in 1:max_depth 
    for dag_vertex in dag.vertices[d]
        operation = dag_vertex.vertex.args[1]
        children = dag_vertex.children 
        children_interval = Any[]
        for child in children
            depth, position = dag_lookup[child]
            push!(children_interval, dag.vertices[depth][position].interval)
        end 
        if length(children_interval) == 1
            new_interval = eval(Expr(:call, operation, children_interval[1]))
            dag_vertex.interval = interserct(dag_vertex.interval, new_interval)
        elseif length(children_interval) == 2 
            new_interval = eval(Expr(:call, operation, children_interval[1], children_interval[2]))
            dag_vertex.interval = intersect(dag_vertex.interval, new_interval)
        else 
            new_expr = Expr(:call, operation, children_interval[1], children_interval[2])
            for i in 3:length(children_interval)
                push!(new_expr.args, children_interval[i])
            end 
            new_interval = eval(new_expr)
            dag_vertex.interval = intersect(dag_vertex.interval, new_interval)
        end 
    end 
end
"""