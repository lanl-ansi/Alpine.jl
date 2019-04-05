__precompile__()

module Alpine

import JuMP
using Memento

using IntervalArithmetic
using IntervalContractors

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

# indices
const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

# functions
const SVF = MOI.SingleVariable
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VECTOR = MOI.VectorOfVariables

# terms 
const SQT = MOI.ScalarQuadraticTerm{Float64}

# sets
const BOUNDS = Union{
    MOI.EqualTo{Float64}, 
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64}, 
    MOI.Interval{Float64}
}

const CONSTRAINT_BOUNDS = Union{ 
    MOI.EqualTo{Float64},
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64}
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

using Compat.Printf 

using LinearAlgebra
using SparseArrays

# Create our module level logger (this will get precompiled)
const LOGGER = getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(Alpine)`
# NOTE: If this line is not included then the precompiled `Alpine.LOGGER` won't be registered at runtime.
__init__() = Memento.register(LOGGER)

"Suppresses information and warning messages output by Alpine, for fine grained control use the Memento package"
function silence()
    info(LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    setlevel!(getlogger(Alpine), "error")
end

include("types.jl")
include("utils.jl")
include("solver_options.jl")
include("logging.jl")
include("expr_type_checking.jl")
include("nl_expr_utilities.jl")
include("fbbt_utilities.jl")
include("fbbt.jl")
include("dag.jl")
include("inits.jl")
include("convexity.jl")
include("user_functions.jl")
include("MOI_wrapper.jl")

end # module
