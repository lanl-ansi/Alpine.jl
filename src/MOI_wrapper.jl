"""
MOI_wrapper.jl defines the Alpine.Optimizer struct
with all mandatory MOI functions overloaded
"""

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
end
function MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x)
    fill!(g, 0.0)
end
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
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
        solver_options)
end 

"""
Printing the optimizer 
"""
function Base.show(io::IO, model::Optimizer)
    println(io, "Alpine Optimizer")
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
        isempty(model.rsoc_constraints)
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
soc_offset(model::Optimizer) = quadratic_eq_offset(model) + length(model.soc_constraints)
rsoc_offset(model::Optimizer) = soc_offset(model) + length(model.rsoc_constraints)
nlp_constraint_offset(model::Optimizer) = rsoc_offset(model) + length(model.quadratic_eq_constraints)

"""
variable offset 
""" 
variable_offset(model::Optimizer) = length(model.variable_info)


"""
``MOI.optimize!()`` for Alpine 
""" 
function MOI.optimize!(model::Optimizer)
    init_ap_data!(model)
    run_convexity_detection!(model)

end 


include(joinpath("MOI_wrapper", "attributes.jl"))
include(joinpath("MOI_wrapper", "variables.jl"))
include(joinpath("MOI_wrapper", "constraints.jl"))
include(joinpath("MOI_wrapper", "objective.jl"))
include(joinpath("MOI_wrapper", "results.jl"))
include(joinpath("MOI_wrapper", "nlp.jl"))