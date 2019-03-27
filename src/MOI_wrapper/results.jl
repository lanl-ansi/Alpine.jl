
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.inner.status.alpine_status == :infeasible 
        return MOI.INFEASIBLE
    end
end