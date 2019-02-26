"""
MOI constraints
"""

"""
Single variable bound constraints 
"""
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.EqualTo{Float64}}) = true

"""
Binary variable support (Integer variables not supported)
"""
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.ZeroOne}) = true

"""
Linear constraints
"""
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{MOI.EqualTo{Float64}}) = true

"""
Quadratic constraints (scalar i.e., vectorized constraints are not supported)
"""
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{MOI.EqualTo{Float64}}) = true

"""
Second-order cone and rotated second-order cone support
"""
MOI.supports_constraint(::Optimizer, ::Type{VECTOR}, ::Type{SOC}) = true
MOI.supports_constraint(::Optimizer, ::Type{VECTOR}, ::Type{RSOC}) = true

"""
``define_add_constraint`` macro - from Ipopt.jl
"""
macro define_add_constraint(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(model::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(model, func)
            push!(model.$(array_name), (func, set))
            return MOI.ConstraintIndex{$function_type, $set_type}(length(model.$(array_name)))
        end
    end
end

"""
``MOI.add_constraint()`` overloads for all the supported constraint types 
"""
@define_add_constraint(SAF, MOI.LessThan{Float64}, linear_le_constraints)
@define_add_constraint(SAF, MOI.GreaterThan{Float64}, linear_ge_constraints)
@define_add_constraint(SAF, MOI.EqualTo{Float64}, linear_eq_constraints)
@define_add_constraint(SQF, MOI.LessThan{Float64}, quadratic_le_constraints)
@define_add_constraint(SQF, MOI.GreaterThan{Float64}, quadratic_ge_constraints)
@define_add_constraint(SQF, MOI.EqualTo{Float64}, quadratic_eq_constraints)
@define_add_constraint(VECTOR, SOC, soc_constraints)
@define_add_constraint(VECTOR, RSOC, rsoc_constraints)

"""
Binary variable support 
"""
function MOI.add_constraint(model::Optimizer, v::SVF, ::MOI.ZeroOne)
    vi = v.variable
	check_inbounds(model, vi)
    model.variable_info[vi.value].lower_bound = 0.0
    model.variable_info[vi.value].has_lower_bound = true
    model.variable_info[vi.value].upper_bound = 1.0
    model.variable_info[vi.value].has_upper_bound = true
    model.variable_info[vi.value].is_binary = true
	
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}(vi.value)
end 

"""
ConstraintIndex support 
"""
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{CI}) = true