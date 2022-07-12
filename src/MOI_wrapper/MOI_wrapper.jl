"""
 MOI_wrapper.jl defines the Alpine.Optimizer struct
 with all mandatory MOI functions overloaded
"""

"""
 Alpine.Optimizer struct
"""
mutable struct Optimizer <: MOI.AbstractOptimizer

    options           :: OptimizerOptions                           # Options set by user

    # Sub-solver identifier for customized solver option
    nlp_solver_id     :: AbstractString                             # NLP Solver identifier string
    minlp_solver_id   :: AbstractString                             # MINLP local solver identifier string
    mip_solver_id     :: AbstractString                             # MIP solver identifier string

    # User inputs
    num_var_orig           :: Int                                   # Initial number of variables
    num_cont_var_orig      :: Int                                   # Initial number of continuous variables
    num_int_var_orig       :: Int                                   # Initial number of binary/integer variables
    num_constr_orig        :: Int                                   # Initial number of constraints
    num_lconstr_orig       :: Int                                   # Initial number of linear constraints
    num_nlconstr_orig      :: Int                                   # Initial number of non-linear constraints
    var_type_orig          :: Vector{Symbol}                        # Variable type vector on original variables (only :Bin, :Cont, :Int)
    var_start_orig         :: Vector{Float64}                       # Variable warm start vector on original variables
    constr_type_orig       :: Vector{Symbol}                        # Constraint type vector on original variables (only :(==), :(>=), :(<=))
    lin_quad_constraints   :: Vector{Any}                           # Constraint `func`-in-`set` values
    constr_expr_orig       :: Vector{Expr}                          # Constraint expressions
    obj_expr_orig          :: Union{Expr,Number}                    # Objective expression

    # Additional user inputs concerning local solves
    l_var_orig                :: Vector{Float64}                    # Variable lower bounds
    u_var_orig                :: Vector{Float64}                    # Variable upper bounds
    constraint_bounds_orig    :: Vector{MOI.NLPBoundsPair}          # Constraint lower bounds
    nl_constraint_bounds_orig :: Vector{MOI.NLPBoundsPair}          # Constraint lower bounds
    sense_orig                :: MOI.OptimizationSense              # Problem type (:Min, :Max)
    d_orig                    :: Union{Nothing, JuMP.NLPEvaluator}  # Instance of AbstractNLPEvaluator for evaluating gradient, Hessian-vector products, and Hessians of the Lagrangian
    disable_hessian           :: Bool                               # Check if there are any user-defined operators, and disable hessians if necessary.
    has_nl_objective          :: Bool
    objective_function        :: Union{Nothing, MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}}

    # Additional initial data
    is_obj_linear_orig        :: Bool                               # Boolean parameter for type of objective

    # Local solution model (extra data for each iteration)
    l_var                     :: Vector{Float64}                    # Updated variable lower bounds for local solve
    u_var                     :: Vector{Float64}                    # Updated variable upper bounds for local solve

    # MIP bounding model
    model_mip                 :: JuMP.Model                         # JuMP's convex-MIP model for bounding

    num_var_linear_mip        :: Int                                # Number of linear lifted variables required
    num_var_nonlinear_mip     :: Int                                # Number of lifted variables
    num_var_disc_mip          :: Int                                # Number of variables which are discretized
    num_constr_convex         :: Int                                # Number of convex constraints

    # Expression-related and structural property placeholder
    linear_terms             :: Dict{Any, Any}                      # Dictionary containing details of lifted linear terms
    nonconvex_terms          :: Dict{Any,Any}                       # Dictionary containing details of lifted non-linear terms
    term_seq                 :: Dict{Int, Any}                      # Vector-Dictionary for NL terms detection
    nonlinear_constrs        :: Dict{Any,Any}                       # Dictionary containing details of special constraints
    obj_structure            :: Symbol                              # A symbolic indicator of the expression type of objective function
    constr_structure         :: Vector{Symbol}                      # A vector indicating whether a constraint is with the special structure
    bounding_obj_expr_mip    :: Union{Expr,Number}                  # Lifted objective expression; if linear, same as obj_expr_orig
    bounding_constr_expr_mip :: Vector{Expr}                        # Lifted constraints; if linear, same as corresponding constr_expr_orig
    bounding_obj_mip         :: Dict{Any, Any}                      # Lifted objective expression in affine form
    bounding_constr_mip      :: Vector{Dict{Any, Any}}              # Lifted constraint expressions in affine form

    # Discretization Related options
    candidate_disc_vars   :: Vector{Int}                            # A vector of all original variable indices that is involved in the nonlinear terms
    discretization        :: Dict{Any,Any}                          # Discretization points with variable keys
    disc_vars             :: Vector{Int}                            # Variables chosen for discretization
    int_vars              :: Vector{Int}                            # Index vector of integer variables
    bin_vars              :: Vector{Int}                            # Index vector of binary variables

    # Reformulated problem options
    l_var_tight           :: Vector{Float64}                        # Post-presolve variable upper bounds
    u_var_tight           :: Vector{Float64}                        # Post-presolve variable lower bounds
    var_type              :: Vector{Symbol}                         # Updated variable type for local solve

    # Solution information
    presolve_infeasible   :: Bool                                   # Presolve infeasibility detection flag
    best_bound            :: Float64                                # Best bound from MIP
    best_obj              :: Float64                                # Best feasible objective value
    initial_warmval       :: Vector{Float64}                        # Warmstart values set to Alpine
    best_sol              :: Vector{Float64}                        # Best feasible solution
    best_bound_sol        :: Vector{Float64}                        # Best bound solution (arg-min)
    presolve_best_rel_gap :: Float64                                # Post-OBBT relative optimality gap = |best_obj - best_bound|/|best_obj|*100
    best_rel_gap          :: Float64                                # Relative optimality gap = |best_obj - best_bound|/|best_obj|*100
    best_abs_gap          :: Float64                                # Absolute gap = |best_obj - best_bound|
    bound_sol_history     :: Vector{Vector{Float64}}                # History of bounding solutions limited by "parameter disc_consecutive_forbid"
    bound_sol_pool        :: Dict{Any, Any}                         # A pool of solutions from solving model_mip

    # Logging information and status
    logs                       :: Dict{Symbol,Any}                          # Logging information
    detected_feasible_solution :: Bool
    detected_bound             :: Bool
    status                     :: Dict{Symbol, MOI.TerminationStatusCode}   # Detailed status of every iteration in the algorithm
    alpine_status              :: MOI.TerminationStatusCode                 # Current Alpine's status

    # Constructor for Alpine.Optimizer
    function Optimizer()
        m = new()
        m.options = Alp.get_default_options()
        MOI.empty!(m)
        return m
    end
end

struct NumberOfIterations <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::NumberOfIterations) = true
MOI.get(m::Optimizer, ::NumberOfIterations)  = m.logs[:n_iter]

struct NumberOfPresolveIterations <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::NumberOfPresolveIterations) = true
MOI.get(m::Optimizer, ::NumberOfPresolveIterations)  = m.logs[:bt_iter]

struct TightenedLowerBound <: MOI.AbstractVariableAttribute end
MOI.is_set_by_optimize(::TightenedLowerBound) = true
function MOI.get(m::Optimizer, ::TightenedLowerBound, vi::MOI.VariableIndex)
    return m.l_var_tight[vi.value]
end

struct TightenedUpperBound <: MOI.AbstractVariableAttribute end
MOI.is_set_by_optimize(::TightenedUpperBound) = true
function MOI.get(m::Optimizer, ::TightenedUpperBound, vi::MOI.VariableIndex)
    return m.u_var_tight[vi.value]
end

MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m.alpine_status

function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    end
    
    status = m.alpine_status
    if status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif status == MOI.NEARLY_FEASIBLE_POINT
        return MOI.NEARLY_FEASIBLE_POINT
    elseif status == MOI.INFEASIBLE_POINT
        return MOI.INFEASIBLE_POINT
    elseif status == MOI.LOCALLY_SOLVED
        return MOI.FEASIBLE_POINT
    end

    return MOI.NO_SOLUTION
end

# MOI.get(::Optimizer, ::MOI.PrimalStatus) = MOI.NO_SOLUTION

MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m.best_obj
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m.best_bound
MOI.get(m::Optimizer, ::MOI.SolveTimeSec)   = m.logs[:total_time]

function get_option(m::Optimizer, s::Symbol)
    getproperty(m.options, s)
end

function set_option(m::Optimizer, s::Symbol, val)
    setproperty!(m.options, s, val)
end

function MOI.is_empty(model::Optimizer)
    return iszero(model.num_var_orig)
end

function MOI.empty!(m::Optimizer)
    m.nlp_solver_id   = ""
    m.minlp_solver_id = ""
    m.mip_solver_id   = ""

    m.num_var_orig         = 0
    m.num_cont_var_orig    = 0
    m.num_int_var_orig     = 0
    m.num_constr_orig      = 0
    m.num_lconstr_orig     = 0
    m.num_nlconstr_orig    = 0
    m.var_type_orig        = Symbol[]
    m.var_start_orig       = Float64[]
    m.constr_type_orig     = Symbol[]
    m.lin_quad_constraints = Any[]
    m.constr_expr_orig     = Expr[]
    # m.num_lconstr_updated = 0
    # m.num_nlconstr_updated = 0
    # m.indices_lconstr_updated = Int[]

    m.l_var_orig                = Float64[]
    m.u_var_orig                = Float64[]
    m.constraint_bounds_orig    = MOI.NLPBoundsPair[]
    m.nl_constraint_bounds_orig = MOI.NLPBoundsPair[]
    m.sense_orig                = MOI.FEASIBILITY_SENSE

    m.d_orig             = nothing
    m.has_nl_objective   = false
    m.objective_function = nothing

    m.linear_terms             = Dict()
    m.nonconvex_terms          = Dict()
    m.term_seq                 = Dict()
    m.nonlinear_constrs        = Dict()
    m.candidate_disc_vars      = Int[]
    m.bounding_constr_expr_mip = []
    m.bounding_constr_mip      = []
    m.disc_vars                = []
    m.int_vars                 = []
    m.bin_vars                 = []
    m.discretization           = Dict()
    m.num_var_linear_mip       = 0
    m.num_var_nonlinear_mip    = 0
    m.num_var_disc_mip         = 0
    m.num_constr_convex        = 0
    m.constr_structure         = Symbol[]
    m.best_bound_sol           = []
    m.bound_sol_history        = []
    m.presolve_infeasible      = false
    m.bound_sol_history        = Vector{Vector{Float64}}(undef, m.options.disc_consecutive_forbid)

    m.best_obj              = Inf
    m.initial_warmval       = Float64[]
    m.best_sol              = Float64[]
    m.best_bound            = -Inf
    m.presolve_best_rel_gap = Inf
    m.best_rel_gap          = Inf
    m.best_abs_gap          = Inf
    m.alpine_status         = MOI.OPTIMIZE_NOT_CALLED

    create_status!(m)
    create_logs!(m)
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOIU.default_copy_to(model, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Alpine"
MOI.get(::Optimizer, ::MOI.SolverVersion) = _ALPINE_VERSION

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    Alp.set_option(model, Symbol(param.name), value)
end
function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    Alp.get_option(model, Symbol(param.name))
end

function MOI.add_variables(model::Optimizer, n::Int)
    return [MOI.add_variable(model) for i in 1:n]
end
function MOI.add_variable(model::Optimizer)
    model.num_var_orig += 1
    push!(model.l_var_orig, -Inf)
    push!(model.u_var_orig, Inf)
    push!(model.var_type_orig, :Cont)
    push!(model.initial_warmval, 0.0)
    push!(model.best_sol, 0.0)
    return MOI.VariableIndex(model.num_var_orig)
end

function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart,
                      ::Type{MOI.VariableIndex})
    return true
end
function MOI.set(model::Optimizer, ::MOI.VariablePrimalStart,
                 vi::MOI.VariableIndex, value::Union{Real, Nothing})
    model.best_sol[vi.value] = model.initial_warmval[vi.value] = something(value, 0.0)
    return
end

const SCALAR_SET = Union{MOI.EqualTo{Float64}, MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.Interval{Float64}}

MOI.supports_constraint(::Optimizer, ::Type{MOI.VariableIndex}, ::Type{<:SCALAR_SET}) = true

_lower(set::MOI.EqualTo)     = set.value
_upper(set::MOI.EqualTo)     = set.value
_lower(set::MOI.LessThan)    = nothing
_upper(set::MOI.LessThan)    = set.upper
_lower(set::MOI.GreaterThan) = set.lower
_upper(set::MOI.GreaterThan) = nothing
_lower(set::MOI.Interval)    = set.lower
_upper(set::MOI.Interval)    = set.upper

function MOI.add_constraint(model::Optimizer, vi::MOI.VariableIndex, set::SCALAR_SET)
    l = _lower(set)
    if l !== nothing
        model.l_var_orig[vi.value] = l
    end
    u = _upper(set)
    if u !== nothing
        model.u_var_orig[vi.value] = u
    end
    return MOI.ConstraintIndex{typeof(vi),typeof(set)}(vi.value)
end

MOI.supports_constraint(::Optimizer, ::Type{MOI.VariableIndex}, ::Type{MOI.Integer}) = true

function MOI.add_constraint(model::Optimizer, f::MOI.VariableIndex, set::MOI.Integer)
    model.var_type_orig[f.value] = :Int
    return MOI.ConstraintIndex{typeof(f), typeof(set)}(f.value)
end
MOI.supports_constraint(::Optimizer, ::Type{MOI.VariableIndex}, ::Type{MOI.ZeroOne}) = true

function MOI.add_constraint(model::Optimizer, f::MOI.VariableIndex, set::MOI.ZeroOne)
    model.var_type_orig[f.value] = :Bin
    return MOI.ConstraintIndex{typeof(f), typeof(set)}(f.value)
end

MOI.supports_constraint(model::Optimizer, ::Type{<:Union{MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}}}, ::Type{<:SCALAR_SET}) = true

function MOI.add_constraint(model::Optimizer, f::Union{MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}}, set::SCALAR_SET)
    model.num_constr_orig += 1
    push!(model.constraint_bounds_orig, MOI.NLPBoundsPair(something(_lower(set), -Inf), something(_upper(set), Inf)))
    iszero(f.constant) || throw(MOI.ScalarFunctionConstantNotZero{Float64, typeof(f), typeof(set)}(f.constant))
    push!(model.lin_quad_constraints, (copy(f), copy(set)))
    push!(model.constr_expr_orig, _constraint_expr(_moi_function_to_expr(f), set))
    if f isa MOI.ScalarAffineFunction
        model.num_lconstr_orig += 1
        push!(model.constr_structure, :generic_linear)
    else
        model.num_nlconstr_orig += 1
        push!(model.constr_structure, :generic_nonlinear)
    end

    return MOI.ConstraintIndex{typeof(f), typeof(set)}(model.num_constr_orig)
end

function MOI.supports(model::Optimizer, ::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction{F}}) where F<:Union{MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}}
    return true
end

is_min_sense(model::Optimizer) = model.sense_orig == MOI.MIN_SENSE
is_max_sense(model::Optimizer) = model.sense_orig == MOI.MAX_SENSE

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense)
    model.sense_orig = sense
    if Alp.is_max_sense(model)
        model.best_obj = -Inf
        model.best_bound = Inf
    elseif Alp.is_min_sense(model)
        model.best_obj = Inf
        model.best_bound = -Inf
    else
        error("Feasibility sense is not supported by Alpine.")
    end
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where F
    model.objective_function = func
end

function MOI.set(m::Optimizer, ::MOI.NLPBlock, block)
    m.d_orig = block.evaluator
    m.disable_hessian = !((:Hess in MOI.features_available(block.evaluator)) && (:HessVec in MOI.features_available(block.evaluator)))
    m.has_nl_objective = block.has_objective

    # We cache it to add it in `load!` as we cannot call `MOI.constraint_expr` yet
    # so we will add the nonlinear `constr_expr_orig` at the end so we need
    # to add the bounds at the end too.
    # So we can consider that the nonlinear constraints are the
    # `length(m.nl_constraint_bounds_orig)` last ones.
    m.nl_constraint_bounds_orig = block.constraint_bounds
end

# In JuMP v0.18/MathProgBase, the 5th decision variable would be `:(x[5])`.
# In JuMP v0.19/MathOptInterface, it is now `:(x[MOI.VariableIndex(5)])`.
# To ease the transition, we simply transform it back to what it used to be.
_variable_index_to_index(expr::Union{Number, Symbol}) = expr
_variable_index_to_index(expr::MOI.VariableIndex) = expr.value
function _variable_index_to_index(expr::Expr)
    for i in eachindex(expr.args)
        expr.args[i] = _variable_index_to_index(expr.args[i])
    end
    return expr
end

_index_to_variable_ref(m::JuMP.Model, idx::Int64) = JuMP.VariableRef(m, MOI.VariableIndex(idx))

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)

    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    
    return model.best_sol[vi.value]
end

MOI.get(model::Optimizer, ::MOI.ResultCount) = model.alpine_status == MOI.OPTIMIZE_NOT_CALLED ? 0 : 1

MOI.is_valid(model::Alpine.Optimizer, vi::MOI.VariableIndex) = 1 <= vi.value <= model.num_var_orig

# Taken from MatrixOptInterface.jl
@enum ConstraintSense EQUAL_TO GREATER_THAN LESS_THAN INTERVAL
_sense(::Type{<:MOI.EqualTo})     = EQUAL_TO
_sense(::Type{<:MOI.LessThan})    = LESS_THAN
_sense(::Type{<:MOI.GreaterThan}) = GREATER_THAN
_sense(::Type{<:MOI.Interval})    = INTERVAL

_no_upper(bound) = bound != typemax(bound)
_no_lower(bound) = bound != typemin(bound)

function _bound_set(lb::T, ub::T) where T
    if _no_upper(ub)
        if _no_lower(lb)
            if ub == lb
                return MOI.EqualTo(lb)
            else
                return MOI.Interval(lb, ub)
            end
        else
            return MOI.LessThan(ub)
        end
    else
        if _no_lower(lb)
            return MOI.GreaterThan(lb)
        else
            return
        end
    end
end

#=
function _bound_sense(lb::T, ub::T) where T
    if _no_upper(ub)
        if _no_lower(lb)
            if lb == ub
                return EQUAL_TO
            else
                return INTERVAL
            end
        else
            return LESS_THAN
        end
    else
        if _no_lower(lb)
            return GREATER_THAN
        else
            return
        end
    end
end
=#