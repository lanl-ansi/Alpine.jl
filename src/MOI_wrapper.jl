using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

""" 
MOI functions, sets and, other type definitions
"""
# indices
const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

# supported functions
const SVF = MOI.SingleVariable
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VECTOR = MOI.VectorOfVariables

# supported sets
const BOUNDS = Union{
    MOI.EqualTo{Float64}, 
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64}, 
    MOI.Interval{Float64}
}
const VAR_TYPES = Union{
    MOI.ZeroOne, 
    MOI.Integer
}
const SOC = MOI.SecondOrderCone
const RSOC = MOI.RotatedSecondOrderCone

# other MOI types
const AFF_TERM = MOI.ScalarAffineTerm{Float64}
const QUAD_TERM = MOI.ScalarQuadraticTerm{Float64}

"""
Variable information struct definition 
""" 
mutable struct VariableInfo
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # Implies lower_bound == Inf
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # Implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound.
    name::String
end
VariableInfo() = VariableInfo(-Inf, false, Inf, false, false, "")

"""
Optimizer struct
"""  
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{AlpineProblem, Nothing}
    variable_info::Vector{VariableInfo}
    nlp_data::MOI.NLPBlockData
    sense::MOI.OptimizationSense 
    objective::Union{SVF, SAF, SQF, Nothing}
    linear_le_constraints::Vector{Tuple{SAF, MOI.LessThan{Float64}}}
    linear_ge_constraints::Vector{Tuple{SAF, MOI.GreaterThan{Float64}}}
    linear_eq_constraints::Vector{Tuple{SAF, MOI.EqualTo{Float64}}}
    quadratic_le_constraints::Vector{Tuple{SAF, MOI.LessThan{Float64}}}
    quadratic_ge_constraints::Vector{Tuple{SQF, MOI.GreaterThan{Float64}}}
    quadratic_eq_constraints::Vector{Tuple{SQF, MOI.EqualTo{Float64}}}
    options
end

# copy constructor
MOIU.supports_default_copy_to(model::Optimizer, copy_names::Bool) = true
function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(model, src; kws...)
end