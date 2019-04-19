"""
    local_solve!(model::MOI.AbstractOptimizer; presolve=false, standalone_solve=false)

The function performs local solve of the original MINLP/NLP model. 
Based on the keyword arguments ``presolve`` and ``standalone_solve``, 
its behaviour will be different. 
"""

function local_solve!(model::MOI.AbstractOptimizer; presolve=false, standalone_solve=false)
    x = getindex(model.inner.continuous_relaxation, :x)
    num_variables = model.inner.num_variables
    binary_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_binary))
    integer_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_integer))
    has_discrete_variables = length(binary_variables) > 0 || length(integer_variables) > 0
    max_multistart_points = model.solver_options.max_multistart_points
    
    if presolve && ~standalone_solve
        if ~isa(get_minlp_optimizer(model), Nothing) && (has_discrete_variables)  
            # set binary and integer variables      
            for i in 1:model.inner.num_variables 
                (binary_variables[i]) && (JuMP.set_binary(x[i]))
                (integer_variables[i]) && (JuMP.set_integer(x[i]))
            end 

            # optimize 
            JuMP.optimize!(model.inner.continuous_relaxation, with_optimizer(get_minlp_optimizer, model))
            status = JuMP.termination_status(model.inner.continuous_relaxation)

            if status_is_optimal(status) 
                model.inner.incumbent.objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                model.inner.incumbent.variable_value = JuMP.value.(x)
                model.inner.incumbent.status = status
                model.inner.status.alpine_status = status
                if model.sense == MOI.MIN_SENSE
                    model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                        -Inf..model.inner.incumbent.objective_value)
                end 
                if model.sense == MOI.MAX_SENSE 
                    model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                        model.inner.incumbent.objective_value..Inf)
                end
            end 

            # unset binary and integer variables
            for i in 1:model.inner.num_variables 
                (binary_variables[i]) && (JuMP.unset_binary(x[i]))
                (integer_variables[i]) && (JuMP.unset_integer(x[i]))
            end  

        elseif isa(get_minlp_optimizer(model), Nothing) && (has_discrete_variables)
            # do nothing 
            info(LOGGER, "Original problem has discrete variables and no minlp_optimizer provided - skipping the local solve in the presolver")
        
        else 
            # optimize the NLP
            JuMP.optimize!(model.inner.continuous_relaxation, with_optimizer(get_nlp_optimizer, model))
            status = JuMP.termination_status(model.inner.continuous_relaxation)

            if status_is_optimal(status) 
                model.inner.incumbent.objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                model.inner.incumbent.variable_value = JuMP.value.(x)
                model.inner.incumbent.status = status
                model.inner.status.alpine_status = status
                if model.sense == MOI.MIN_SENSE
                    model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                        -Inf..model.inner.incumbent.objective_value)
                end 
                if model.sense == MOI.MAX_SENSE 
                    model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                        model.inner.incumbent.objective_value..Inf)
                end
            end 

        
            # one random restart for alternate local optimal solution 
            start_values = generate_random_start_values(model)
            for k in 1:num_variables 
                JuMP.set_start_value(x[k], start_values[k])
            end 
            JuMP.optimize!(model.inner.continuous_relaxation)
            status = JuMP.termination_status(model.inner.continuous_relaxation)

            if status_is_optimal(status) 
                current_objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                is_better = false 
                (model.sense == MOI.MIN_SENSE) && (is_better = (current_objective_value < model.inner.incumbent.objective_value - 1e-4))
                (model.sense == MOI.MAX_SENSE) && (is_better = (current_objective_value > model.inner.incumbent.objective_value + 1e-4))
                if is_better 
                    model.inner.incumbent.objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                    model.inner.incumbent.variable_value = JuMP.value.(x)
                    model.inner.incumbent.status = status
                    model.inner.status.alpine_status = status
                    if model.sense == MOI.MIN_SENSE
                        model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                            -Inf..model.inner.incumbent.objective_value)
                    end 
                    if model.sense == MOI.MAX_SENSE 
                        model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                            model.inner.incumbent.objective_value..Inf)
                    end
                end
            end
            
        end 
    
    elseif presolve && standalone_solve && ~has_discrete_variables
        for i in 1:max_multistart_points
            start_values = generate_random_start_values(model)
            for k in 1:num_variables 
                JuMP.set_start_value(x[k], start_values[k])
            end 
            JuMP.optimize!(model.inner.continuous_relaxation)

            status = JuMP.termination_status(model.inner.continuous_relaxation)

            if status_is_optimal(status) 
                current_objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                is_better = false 
                (model.sense == MOI.MIN_SENSE) && (is_better = (current_objective_value < model.inner.incumbent.objective_value - 1e-4))
                (model.sense == MOI.MAX_SENSE) && (is_better = (current_objective_value > model.inner.incumbent.objective_value + 1e-4))
                if is_better 
                    model.inner.incumbent.objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                    model.inner.incumbent.variable_value = JuMP.value.(x)
                    model.inner.incumbent.status = status
                    model.inner.status.alpine_status = status
                    if model.sense == MOI.MIN_SENSE
                        model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                            -Inf..model.inner.incumbent.objective_value)
                    end 
                    if model.sense == MOI.MAX_SENSE 
                        model.inner.objective_bound_info = intersect(model.inner.objective_bound_info, 
                            model.inner.incumbent.objective_value..Inf)
                    end
                end
            end
        end

    else 
        # called during bounding solve iterations (complete later)
    
    end 
    return

end