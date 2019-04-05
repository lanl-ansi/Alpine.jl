
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if isa(model.inner.status, Nothing)
        return MOI.NO_SOLUTION
    elseif model.inner.status.alpine_status == :infeasible 
        return MOI.INFEASIBLE
    end
end
