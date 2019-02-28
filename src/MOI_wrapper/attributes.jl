"""
``is_convex`` constraint attribute - this is holds only for quadratic constraints
"""

struct is_convex <: MOI.AbstractConstraintAttribute end

MOI.supports(::Optimizer, is_convex::MOI.AbstractConstraintAttribute, ::Type{CI{SQF, MOI.LessThan{Float64}}}) = true
MOI.supports(::Optimizer, is_convex::MOI.AbstractConstraintAttribute, ::Type{CI{SQF, MOI.GreaterThan{Float64}}}) = true

MOI.is_copyable(::is_convex) = true