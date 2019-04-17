"""
    local_solve!(model::MOI.AbstractOptimizer; presolve=false)

The function performs local solve of the original MINLP/NLP model. 
If the original problem is an MINLP and if an MINLP solver is 
provided using the ``minlp_optimizer`` option, then an MINLP local
solve is performed at the presolve step; if no MINLP solver is 
provided, then no local solve is performed in the presolve step 
(this would impact the bound-tightening algorithm and result in 
weaker variable bounds). If the original problem is an NLP, the 
presolve step always performs a local solve with the NLP solver 
(mandatory option) that was provided using the ``nlp_optimzer`` 
option. In the local solves performed during each iteration of 
the bounding solve, Alpine.jl always solves an NLP to local 
optimality. When the original problem is an MINLP, the discrete 
variables get fixed to the value obtained at that corresponding 
iteration of the bounding solve. 
"""

function local_solve!(model::MOI.AbstractOptimizer; presolve=false)
    x = getindex(model.inner.continuous_relaxation, :x)
    num_variables = model.inner.num_variables
    binary_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_binary))
    integer_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_integer))
    has_discrete_variables = length(binary_variables) > 0 || length(integer_variables) > 0
    max_multistart_points = model.solver_options.max_multistart_points * 0

    if presolve 
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
            @show MOI.is_empty(model.solver_options.nlp_optimizer)
            JuMP.optimize!(model.inner.continuous_relaxation, with_optimizer(get_nlp_optimizer, model))
            status = JuMP.termination_status(model.inner.continuous_relaxation)

            if status_is_optimal(status) 
                model.inner.incumbent.objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                model.inner.incumbent.variable_value = JuMP.value.(x)
                model.inner.incumbent.status = status
                model.inner.status.alpine_status = status
            end 

            """
        
            # multi-start NLP for alternate local optimal solutions 
            for i in 1:max_multistart_points
                start_values = generate_random_start_values(model)
                for k in 1:num_variables 
                    JuMP.set_start_value(x[k], start_values[k])
                end 
                MOIU.drop_optimizer(JuMP.backend(model.inner.continuous_relaxation))
                JuMP.optimize!(model.inner.continuous_relaxation, with_optimizer(get_nlp_optimizer, model))
                status = JuMP.termination_status(model.inner.continuous_relaxation)

                if status_is_optimal(status) 
                    current_objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                    is_better = false 
                    (model.sense == MOI.MIN_SENSE) && (is_better = (current_objective_value < model.inner.incumbent.objective_value - 1e-4))
                    (model.sense == MOI.MAX_SENSE) && (is_vetter = (current_objective_value > model.inner.incumbent.objective_value + 1e-4))
                    if is_better 
                        model.inner.incumbent.objective_value = JuMP.objective_value(model.inner.continuous_relaxation) 
                        model.inner.incumbent.variable_value = JuMP.value.(x)
                        model.inner.incumbent.status = status
                        model.inner.status.alpine_status = status
                    end
                end
            end

            """  
        end 
    
    else 
        # called during bounding solve iterations (complete later)
    
    end 
    return

end