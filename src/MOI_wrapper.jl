"""
MOI_wrapper.jl defines the Alpine.Optimizer struct
with all mandatory MOI functions overloaded
"""

""" 
MOI functions, sets and, other type definitions
"""
# indices
const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

# functions
const SVF = MOI.SingleVariable
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VECTOR = MOI.VectorOfVariables

# sets
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
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound
    is_binary::Bool       # Implies lower_bound == 0, upper_bound == 1 and is MOI.ZeroOne
    name::String
end
VariableInfo() = VariableInfo(-Inf, false, Inf, false, false, false, "")

"""
function to get an array of variable attributes 
"""
function info_array_of_variables(variable_info::Vector{VariableInfo}, attr::Symbol)
    len_var_info = length(variable_info)
    type_dict = get_type_dict(variable_info[1])
    result = Array{type_dict[attr], 1}(undef, len_var_info)
    for i = 1:len_var_info
        result[i] = getfield(variable_info[i], attr)
    end
    return result
end

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
    quadratic_le_constraints::Vector{Tuple{SQF, MOI.LessThan{Float64}}}
    quadratic_ge_constraints::Vector{Tuple{SQF, MOI.GreaterThan{Float64}}}
    quadratic_eq_constraints::Vector{Tuple{SQF, MOI.EqualTo{Float64}}}
    soc_constraints::Vector{Tuple{VECTOR, SOC}}
    rsoc_constraints::Vector{Tuple{VECTOR, RSOC}}
    convex_constraint_indices::Vector{CI}
    solver_options::SolverOptions
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Alpine"

"""
EmptyNLPEvaluator struct and associated functions 
"""
struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end
MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
function MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x)
    fill!(g, 0.0)
    return
end
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end
empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

"""
Optimizer struct constructor 
"""
function Optimizer(; options...) 
    solver_options = combine_options(options)

    return Optimizer(
        nothing, 
        [], 
        empty_nlp_data(), 
        MOI.FEASIBILITY_SENSE, 
        nothing, 
        [], [], [], # linear constraints 
        [], [], [], # quadratic constraints
        [], [], # conic constraints 
        [], # convex constraint ids
        solver_options)
end 

function Optimizer(options::Dict{Symbol,Any}) 
    solver_options = combine_options(options)

    return Optimizer(
        nothing, 
        [], 
        empty_nlp_data(), 
        MOI.FEASIBILITY_SENSE, 
        nothing, 
        [], [], [], # linear constraints 
        [], [], [], # quadratic constraints
        [], [], # conic constraints 
        [], # convex constraint ids
        solver_options)
end 

"""
Printing the optimizer 
"""
function Base.show(io::IO, model::Optimizer)
    println(io, "Alpine Optimizer")
    return
end

""" 
Copy constructor for the optimizer
"""
MOIU.supports_default_copy_to(model::Optimizer, copy_names::Bool) = true
function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(model, src; kws...)
end

"""
``MOI.is_empty(model::Optimizer)`` overload for Alpine.Optimizer
"""
function MOI.is_empty(model::Optimizer)
    return isempty(model.variable_info) &&
        model.nlp_data.evaluator isa EmptyNLPEvaluator &&
        model.sense == MOI.FEASIBILITY_SENSE &&
        isempty(model.linear_le_constraints) &&
        isempty(model.linear_ge_constraints) &&
        isempty(model.linear_eq_constraints) &&
        isempty(model.quadratic_le_constraints) &&
        isempty(model.quadratic_ge_constraints) &&
        isempty(model.quadratic_eq_constraints) && 
        isempty(model.soc_constraints) && 
        isempty(model.rsoc_constraints) &&
        isempty(model.convex_constraint_indices)
end

"""
``MOI.empty!(model::Optimizer)`` overload for Alpine.Optimizer
"""
function MOI.empty!(model::Optimizer)
    model.inner = nothing
    empty!(model.variable_info)
    model.nlp_data = empty_nlp_data()
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
    empty!(model.linear_le_constraints)
    empty!(model.linear_ge_constraints)
    empty!(model.linear_eq_constraints)
    empty!(model.quadratic_le_constraints)
    empty!(model.quadratic_ge_constraints)
    empty!(model.quadratic_eq_constraints)
    empty!(model.soc_constraints)
    empty!(model.rsoc_constraints)
    empty!(model.convex_constraint_indices)
end

"""
ordering of constraints provided to Alpine.jl 
"""
linear_le_offset(model::Optimizer) = 0
linear_ge_offset(model::Optimizer) = length(model.linear_le_constraints)
linear_eq_offset(model::Optimizer) = linear_ge_offset(model) + length(model.linear_ge_constraints)
quadratic_le_offset(model::Optimizer) = linear_eq_offset(model) + length(model.linear_eq_constraints)
quadratic_ge_offset(model::Optimizer) = quadratic_le_offset(model) + length(model.quadratic_le_constraints)
quadratic_eq_offset(model::Optimizer) = quadratic_ge_offset(model) + length(model.quadratic_ge_constraints)
nlp_constraint_offset(model::Optimizer) = quadratic_eq_offset(model) + length(model.quadratic_eq_constraints)

"""
populate_generic_ap_data!() for Alpine 
""" 
function populate_generic_ap_data!(model::Optimizer)
    model.inner = AlpineProblem()
    model.inner.num_variables = length(model.variable_info)
    model.inner.num_linear_eq_constraints = length(model.linear_le_constraints)
    model.inner.num_linear_ge_constraints = length(model.linear_ge_constraints)
    model.inner.num_linear_eq_constraints = length(model.linear_eq_constraints)
    model.inner.num_quadratic_le_constraints = length(model.quadratic_le_constraints)
    model.inner.num_quadratic_ge_constraints = length(model.quadratic_ge_constraints)
    model.inner.num_quadratic_eq_constraints = length(model.quadratic_eq_constraints)
    model.inner.num_soc_constraints = length(model.soc_constraints)
    model.inner.num_rsoc_constraints = length(model.rsoc_constraints) 

    model.inner.lower_original = Vector{Float64}()
    model.inner.upper_original = Vector{Float64}()

    for vi in model.variable_info
        if (!vi.has_lower_bound || !vi.has_upper_bound) 
            error(LOGGER, "Alpine.jl requires every variable in the problem to be bounded; 
                ensure that bounds are provided for every variable")
        end 
        push!(model.inner.lower_original, vi.lower_bound)
        push!(model.inner.upper_original, vi.upper_bound)
    end 

    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)
        num_nlp_constraints = length(model.nlp_data.constraint_bounds)
        model.inner.num_nlp_constraints = num_nlp_constraints
    else 
        info(LOGGER, "no explicit NLP constraints or objective provided using @NLconstraint or @NLobjective macros")
        model.inner.num_nlp_constraints = 0
    end 

    model.inner.num_constraints = model.inner.num_linear_eq_constraints + 
        model.inner.num_linear_ge_constraints +
        model.inner.num_linear_eq_constraints +
        model.inner.num_quadratic_le_constraints +
        model.inner.num_quadratic_ge_constraints +
        model.inner.num_quadratic_eq_constraints +
        model.inner.num_soc_constraints +
        model.inner.num_rsoc_constraints + 
        model.inner.num_nlp_constraints
    
    print_var_con_summary(model.inner)
end

"""
``MOI.optimize!()`` for Alpine 
""" 
function MOI.optimize!(model::Optimizer)
    populate_generic_ap_data!(model)

end 


include(joinpath("MOI_wrapper", "variables.jl"))
include(joinpath("MOI_wrapper", "constraints.jl"))
include(joinpath("MOI_wrapper", "objective.jl"))
include(joinpath("MOI_wrapper", "results.jl"))
include(joinpath("MOI_wrapper", "nlp.jl"))