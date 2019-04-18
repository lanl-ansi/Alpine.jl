
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.inner.status.alpine_status == MOI.OPTIMIZE_NOT_CALLED
        return MOI.NO_SOLUTION 
    end
    return model.inner.status.alpine_status
end
