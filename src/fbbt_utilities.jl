"""
Forward propagation for linear constraints 
""" 
function fbbt_forward_linear_constraints!(model::MOI.AbstractOptimizer)

    for i in 1:model.inner.num_linear_constraints 
        func, set = model.linear_constraints[i]
        offset = linear_offset(model)
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

    return

end 


"""
Backward propagation for linear constraints 
"""
function fbbt_backward_linear_constraints!(model::MOI.AbstractOptimizer)

    # Reference: https://doi.org/10.1080/10556788.2017.1335312 (section 2.1)
    for i in 1:model.inner.num_linear_constraints
        func, set = model.linear_constraints[i]
        offset = linear_offset(model)
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

    for i in 1:model.inner.num_quadratic_constraints 
        func, set = model.quadratic_constraints[i]
        offset = quadratic_offset(model)
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
                    computed_bound_info += coeff..coeff 
                elseif fname == :linear_part 
                    var_id = expr.args[2].value
                    computed_bound_info += coeff * model.inner.variable_bound_tightened[var_id]
                else 
                    depth, position = depth, position = dag_lookup[expr]
                    computed_bound_info += coeff * dag.vertices[depth][position].interval
                end 
                @show expr, computed_bound_info
            end 
        end 
        @show model.inner.constraint_bound_info[offset + i], computed_bound_info
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
                push!(children_interval, dag.vertices[depth][position].interval)
            end 
            propagate_intervals!(operation, current_interval, children_interval)
            i = 1
            for child in children 
                depth, position = dag_lookup[child]
                dag.vertices[depth][position].interval = children_interval[i]
                @show dag.vertices[depth][position].interval
                i += 1
            end
        end 
    end 


    return
    
end 

"""
Propagate intervals backward given an operation, current_interval and intervals for each child 
Reference: https://edoc.hu-berlin.de/bitstream/handle/18452/17356/vigerske.pdf?sequence=1&isAllowed=y (Page 175)
"""
function propagate_intervals!(operation::Symbol, current_interval::Interval{Float64}, children_interval::Vector{Interval{Float64}})

    # a = -b 
    if operation == :- && length(children_interval) == 1
        b = children_interval[1] ∩ (-current_interval)
        children_interval[1] = b
    # a = b - c
    elseif operation == :- && length(children_interval) == 2 
        b = children_interval[1] ∩ (current_interval + children_interval[2])
        c = children_interval[2] ∩ (children_interval[1] - current_interval)
        children_interval[1] = b 
        children_interval[2] = c
    # a = b + c
    elseif operation == :+ && length(children_interval) == 2
        b = children_interval[1] ∩ (current_interval - children_interval[2])
        c = children_interval[2] ∩ (current_interval - children_interval[1])
        children_interval[1] = b 
        children_interval[2] = c
    # a = b * c
    elseif operation == :* && length(children_interval) == 2
        ((0.0 ∉ current_interval) || (0.0 ∉ children_interval[1])) && (children_interval[2] = children_interval[2] ∩ (current_interval / children_interval[1]))
        ((0.0 ∉ current_interval) || (0.0 ∉ children_interval[2])) && (children_interval[1] = children_interval[1] ∩ (current_interval / children_interval[2]))
    # a = b / c
    elseif operation == :/ 
        children_interval[1] = children_interval[1] ∩ (current_interval * children_interval[2])
        children_interval[2] = children_interval[2] ∩ (children_interval[1] / current_interval)
    # a = exp(b)
    elseif operation == :exp 
        lo = max(current_interval.lo, 0.0)
        hi = current_interval.hi
        children_interval[1] = children_interval[1] ∩ log(lo..hi) 
    # a = log(b)
    elseif operation == :log 
        children_interval[1] = children_interval[1] ∩ exp(current_interval) 
    # a = sqrt(b)
    elseif operation == :sqrt 
        lo = max(current_interval.lo, 0.0)
        hi = current_interval.hi
        children_interval[1] = children_interval[1] ∩ (lo..hi)^2
    # a = b^c (assuming c is always a constant interval) -- implement c::Interval{Float64} later
    elseif operation == :^ && isthin(children_interval[2])
        a = current_interval 
        b = children_interval[1]
        c = children_interval[2].lo
        is_integer = try isa(Int(c), Int) catch y false end 
        if is_integer 
            if n == 2  
                root = sqrt(a)
                b1 = b ∩ root
                b2 = b ∩ (-root)
            elseif iseven(n)
                root = a^(1//n)
                b1 = b ∩ root
                b2 = b ∩ (-root)
            elseif isodd(n)
                pos_root = (a ∩ (0..Inf)) ^ (1//n)
                neg_root = -( ( (-a) ∩ (0..Inf) ) ^ (1//n) )
                b1 = b ∩ pos_root
                b2 = b ∩ neg_root
            end
            children_interval[1] = hull(b1, b2)
        else
            lo = max(a.lo, 0) 
            hi = a.hi 
            children_interval[1] = children_interval[1] ∩ lo..hi^(1//c)
            info(LOGGER, "Lower bound of a variable/expression with fractional power set to 0.0 to avoid complex number results")
        end
    # a = abs(b)
    elseif operation == :abs 
        a = current_interval 
        b = children_interval[1]
        a_new = a ∩ (0..Inf)
        b1 = a_new ∩ b
        b2 = -(a_new ∩ (-b))
        children_interval[1] = hull(x1, x2)
    # trigonometric functions
    elseif operation == :sin 
        current_interval, children_interval[1] = sin_rev(current_interval, children_interval[1])
    elseif operation == :cos 
        current_interval, children_interval[1] = cos_rev(current_interval, children_interval[1])
    elseif operation == :tan 
        current_interval, children_interval[1] = tan_rev(current_interval, children_interval[1])
    # k-nary addition 
    elseif operation == :+ && length(children_interval) > 2 
        for i in 1:length(children_interval)
            sum_interval = 0..0
            for j in 1:length(children_interval)
                (j == i) && (continue)
                sum_interval += children_interval[j]
            end 
            children_interval[i] = children_interval[i] ∩ (current_interval - sum_interval)
        end 
    # k-nary multiplication
    elseif operation == :* && length(children_interval) > 2 
        for i in 1:length(children_interval)
            product_interval = 1..1
            for j in 1:length(children_interval)
                (j == i) && (continue) 
                product_interval *= children_interval[j]
            end 
            ((0.0 ∉ current_interval) || (0.0 ∉ product_interval)) && (children_interval[i] = children_interval[i] ∩ (current_interval / product_interval))
        end 
    end 

    return 
end 

"""
Backward propogation for quadratic constraints
"""
function fbbt_backward_quadratic_constraints!(model::MOI.AbstractOptimizer)
    fnames = [:constant_part, :linear_part, :quadratic_part, :bilinear_part]
    dag_lookup = model.inner.dag_lookup
    dag = model.inner.expression_graph
    for i in 1:model.inner.num_quadratic_constraints 
        nl_function = model.inner.quadratic_nl_function[i]
        offset = quadratic_offset(model)
        constraint_interval = model.inner.constraint_bound_info[offset + i] 
        for fname in fnames
            (fname == :constant_part) && (continue)  
            f_value = getfield(nl_function, fname)
            (isa(f_value, Nothing)) && (continue)
            for alpine_expression in f_value
                coeff, expr = alpine_expression.expression
                depth, position = dag_lookup[expr]
                sum_interval = 0..0
                for fname_1 in fnames 
                    f_value_1 = getfield(nl_function, fname_1)
                    (isa(f_value_1, Nothing)) && (continue)
                    for sum_expression in f_value_1 
                        coeff_1, expr_1 = sum_expression.expression
                        (fname == :constant_part) && (sum_interval += coeff_1..coeff_1; continue)
                        (~haskey(dag_lookup, expr_1)) && (error(LOGGER, "expression $expr_1 not found in expression graph"))
                        depth_1, position_1 = dag_lookup[expr_1]
                        if depth_1 == depth && position_1 == position
                            continue
                        end 
                        sum_interval = sum_interval + coeff_1..coeff * dag.vertices[depth][position].interval
                    end 
                end 
                dag.vertices[depth][position].interval = dag.vertices[depth][position].interval ∩ (1/coeff * (constraint_interval - sum_interval))
            end
        end 
    end 

    return
end 

"""
Backward propogation for nonlinear constraints
"""
function fbbt_backward_nl_constraints!(model::MOI.AbstractOptimizer)
    fnames = fieldnames(NLFunction)
    dag_lookup = model.inner.dag_lookup
    dag = model.inner.expression_graph
    for i in 1:model.inner.num_quadratic_constraints 
        nl_function = model.inner.quadratic_nl_function[i]
        offset = quadratic_offset(model)
        constraint_interval = model.inner.constraint_bound_info[offset + i] 
        for fname in fnames
            (fname == :constant_part) && (continue)  
            f_value = getfield(nl_function, fname)
            (isa(f_value, Nothing)) && (continue)
            for alpine_expression in f_value
                coeff, expr = alpine_expression.expression
                depth, position = dag_lookup[expr]
                sum_interval = 0..0
                for fname_1 in fnames 
                    f_value_1 = getfield(nl_function, fname_1)
                    (isa(f_value_1, Nothing)) && (continue)
                    for sum_expression in f_value_1 
                        coeff_1, expr_1 = sum_expression.expression
                        (fname == :constant_part) && (sum_interval += coeff_1..coeff_1; continue)
                        (~haskey(dag_lookup, expr_1)) && (error(LOGGER, "expression $expr_1 not found in expression graph"))
                        depth_1, position_1 = dag_lookup[expr_1]
                        if depth_1 == depth && position_1 == position
                            continue
                        end 
                        sum_interval = sum_interval + coeff_1..coeff * dag.vertices[depth][position].interval
                    end 
                end 
                dag.vertices[depth][position].interval = dag.vertices[depth][position].interval ∩ (1/coeff * (constraint_interval - sum_interval))
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

