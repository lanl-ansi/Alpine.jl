
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    return model.inner.status.alpine_status
end
