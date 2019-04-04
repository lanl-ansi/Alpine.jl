
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.inner.status.alpine_status == :infeasible 
        return MOI.INFEASIBLE
    end
end

function lower_bound_tightened(vref::JuMP.VariableRef)
    index = vref.index.value 
    return vref.model.moi_backend.optimizer.model.inner.variable_bound_tightened[index].lo
end 

function upper_bound_tightened(vref::JuMP.VariableRef)
    index = vref.index.value 
    return vref.model.moi_backend.optimizer.model.inner.variable_bound_tightened[index].hi
end 