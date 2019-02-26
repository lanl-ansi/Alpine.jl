"""
MOI objective
"""

"""
Supported objective types and sense 
"""
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SVF}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SAF}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SQF}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
    
"""
set and get function overloads
"""
MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.sense = sense
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction, func::Union{SVF, SAF, SQF})
    check_inbounds(model, func)
    model.objective = func
end