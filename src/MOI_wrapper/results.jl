
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.inner.status.alpine_status == MOI.OPTIMIZE_NOT_CALLED
        return MOI.NO_SOLUTION 
    end
    return model.inner.status.alpine_status
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    if model.inner.status.alpine_status == MOI.OPTIMIZE_NOT_CALLED 
        error(LOGGER, "Either the JuMP.Optimize! function has not been called on the model or it has been called with one of the options `perform_bounding_solve_only`, `perform_bp_only` set to true")
    end
    return model.inner.incumbent.objective_value
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if model.inner.status.alpine_status == MOI.OPTIMIZE_NOT_CALLED
        error(LOGGER, "Either the JuMP.Optimize! function has not been called on the model or it has been called with one of the options `perform_bounding_solve_only`, `perform_bp_only` set to true")
    end
    check_inbounds(model, vi)
    return model.inner.incumbent.variable_value[vi.value]
end

