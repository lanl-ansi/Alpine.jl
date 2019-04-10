"""
Functions to build the continuous relaxation of the model if the model is an MINLP; 
if not we build the NLP for local solving 
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

    # add all the constraints (convert SOC and RSOC constraints using bridges)

    
    return 
end