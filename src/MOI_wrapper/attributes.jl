"""
``is_convex`` constraint attribute - this is holds only for quadratic constraints
"""

struct is_convex <: MOI.AbstractConstraintAttribute end

MOI.supports(::Optimizer, is_convex::MOI.AbstractConstraintAttribute, ::Type{CI{SQF, MOI.LessThan{Float64}}}) = true
MOI.supports(::Optimizer, is_convex::MOI.AbstractConstraintAttribute, ::Type{CI{SQF, MOI.GreaterThan{Float64}}}) = true

MOI.is_copyable(::is_convex) = true

"""
Variable attributes for getting tightened lower and upper bounds after bound-tightening step 
"""

struct VariableLBTightened <: MOI.AbstractVariableAttribute end 
struct VariableUBTightened <: MOI.AbstractVariableAttribute end 
struct VariableBoundTightened <: MOI.AbstractVariableAttribute end 

MOI.is_set_by_optimize(::Union{VariableLBTightened, VariableUBTightened, VariableBoundTightened}) = true 
MOI.supports(::Optimizer, ::Union{VariableLBTightened, VariableUBTightened, VariableBoundTightened}, ::Type{VI}) = true

function MOI.get(model::Optimizer, attr::VariableLBTightened, vi::VI) 
    return model.inner.variable_bound_tightened[vi.value].lo
end 

function MOI.get(model::Optimizer, attr::VariableUBTightened, vi::VI) 
    return model.inner.variable_bound_tightened[vi.value].hi
end 

function MOI.get(model::Optimizer, attr::VariableBoundTightened, vi::VI)
    interval = model.inner.variable_bound_tightened[vi.value]
    return interval
end

function lower_bound_tightened(v::JuMP.VariableRef)::Float64
    return MOI.get(JuMP.owner_model(v), VariableLBTightened(), v)
end

function upper_bound_tightened(v::JuMP.VariableRef)::Float64
    return MOI.get(JuMP.owner_model(v), VariableUBTightened(), v)
end

 
function variable_bound_tightened(v::JuMP.VariableRef)::Interval{Float64}
    return MOI.get(JuMP.owner_model(v), VariableBoundTightened(), v)
end

"""
Model attributes 
"""
struct IsProblemConvex <: MOI.AbstractModelAttribute end

MOI.is_set_by_optimize(::IsProblemConvex) = true 
MOI.supports(::Optimizer, ::IsProblemConvex) = true 
MOI.get(model::Optimizer, attr::IsProblemConvex) = model.inner.is_problem_convex 
is_problem_convex(m::JuMP.Model)::Bool = MOI.get(m, IsProblemConvex())
