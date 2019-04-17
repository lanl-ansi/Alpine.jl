"""
    build_local_solve_model!(model::MOI.AbstractOptimizer)

The function builds a NLP to provide it to the local solver or in the case 
of an MINLP, it builds a continuous relaxation. 
"""

function build_local_solve_model!(model::MOI.AbstractOptimizer)
    m = JuMP.Model()

    num_variables = model.inner.num_variables
    warm_start = info_array_of_variables(model.variable_info, :start)
    lb = [ bound.lo for bound in model.inner.variable_bound_tightened ]
    ub = [ bound.hi for bound in model.inner.variable_bound_tightened ]

    @variable(m, lb[i] <= x[i=1:num_variables] <= ub[i])
    
    # set start values
    for i in 1:num_variables 
        if ~isa(warm_start[i], Nothing)
            JuMP.set_start_value(x[i], warm_start[i])
        end 
    end 

    backend = JuMP.backend(m)

    # add linear constraints 
    linear_constraints = model.linear_constraints 
    for constraint in linear_constraints 
        MOI.add_constraint(backend, constraint[1], constraint[2])
    end 

    # add quadratic constraints 
    quadratic_constraints = model.quadratic_constraints 
    for constraint in quadratic_constraints 
        MOI.add_constraint(backend, constraint[1], constraint[2])
    end 

    # add SOC and RSOC constraints 
    soc_constraints = model.soc_constraints 
    rsoc_constraints = model.rsoc_constraints

    # add nonlinear constraints 
    for i in 1:model.inner.num_nlp_constraints
        constraint_expr = MOI.constraint_expr(model.nlp_data.evaluator, i)
        dereference_expr!(constraint_expr, m)
        JuMP.add_NL_constraint(m, constraint_expr)
    end

    # add objective 
    if model.inner.is_objective_nl 
        objective_expr = MOI.objective_expr(model.nlp_data.evaluator)
        dereference_expr!(objective_expr, m)
        JuMP.set_NL_objective(m, model.sense, objective_expr)
    elseif ~isa(model.objective, Nothing)
        MOI.set(m, MOI.ObjectiveFunction{typeof(model.objective)}(), model.objective)
        MOI.set(m, MOI.ObjectiveSense(), model.sense)
    end

    model.inner.continuous_relaxation = m
    
    return 
end

""" 
    generate_random_start_values(model::AbstractOptimizer)::Vector{Float64}

This function generates a random warm start value for each variable in the problem. 
It returns a vector of warm start values that is used later in a multi-start setting.
"""
function generate_random_start_values(model::MOI.AbstractOptimizer)::Vector{Float64}
    num_variables = model.inner.num_variables 
    start_values = Vector{Float64}(undef, num_variables)

    for i in 1:num_variables 
        lb = model.inner.variable_bound_tightened[i].lo 
        ub = model.inner.variable_bound_tightened[i].hi 

        start_values[i] = round(rand() * (ub - lb) + lb, digits = 2)
    end 

    return start_values
end 