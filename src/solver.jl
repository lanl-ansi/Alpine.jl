mutable struct OptimizerOptions
    # Parameters for tuning Alpine

    # Basic solver parameters
    log_level::Int                                              # Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info
    time_limit::Float64                                         # Time limit for algorithm (in seconds)
    max_iter::Int                                               # Target Maximum Iterations
    rel_gap::Float64                                            # Relative optimality gap for termination condition
    gap_ref::Symbol                                             # Relative gap reference point (options: [:ub, :lb])
    abs_gap::Float64                                            # Absolute optimality gap for termination condition
    tol::Float64                                                # Numerical tol used in algorithms
    large_bound::Float64                                        # Large bounds for problems with unbounded variables

    # All the solver options
    nlp_solver                                                  # Local continuous NLP solver for solving NLPs at each iteration
    minlp_solver                                                # Local MINLP solver for solving MINLPs at each iteration
    mip_solver                                                  # MIP solver for successive lower bound solves

    # Convexification methods
    recognize_convex::Bool                                      # Recognize convex expressions in parsing objective functions and constraints
    bilinear_mccormick::Bool                                    # [INACTIVE] Convexify bilinear terms using piecwise McCormick representation
    bilinear_convexhull::Bool                                   # Convexify bilinear terms using lambda representation
    monomial_convexhull::Bool                                   # Convexify monomial terms using convex-hull representation

    # Expression-based user-inputs
    method_convexification::Array{Function}                     # Array of functions that user can choose to convexify specific non-linear terms : no over-ride privilege

    # method_partition_injection::Array{Function}               # [INACTIVE] Array of functions for special methods to add partitions to variables under complex conditions
    # term_patterns::Array{Function}                            # [INACTIVE] Array of functions that user can choose to parse/recognize nonlinear terms in constraint expressions
    # constr_patterns::Array{Function}                          # [INACTIVE] Array of functions that user can choose to parse/recognize structural constraint from expressions

    # Parameters used in Alpine's MIP-based partitioning algorithm
    disc_var_pick::Any                                          # Algorithm for choosing the variables to discretize: 1 for minimum vertex cover, 0 for all variables
    disc_ratio::Any                                             # Discretization ratio parameter (use a fixed value for now, later switch to a function)
    disc_uniform_rate::Int                                      # Discretization rate parameter when using uniform partitions
    disc_add_partition_method::Any                              # Additional methods to add discretization
    disc_divert_chunks::Int                                     # How many uniform partitions to construct
    disc_abs_width_tol::Float64                                 # Absolute tolerance used when setting up partition/discretization
    disc_rel_width_tol::Float64                                 # Relative width tolerance when setting up partition/discretization
    disc_consecutive_forbid::Int                                # Prevent bounding model to add partitions consecutively in the same region when bounds do not improve
    disc_ratio_branch::Bool                                     # Branching tests for picking fixed discretization ratio

    # MIP formulation parameters
    convhull_formulation::String                                # MIP Formulation for the relaxation
    convhull_ebd::Bool                                          # Enable embedding formulation
    convhull_ebd_encode::Any                                    # Encoding method used for convhull_ebd
    convhull_ebd_ibs::Bool                                      # Enable independent branching scheme
    convhull_ebd_link::Bool                                     # Linking constraints between x and α, type 1 uses hierarchical and type 2 uses big-m
    convhull_warmstart::Bool                                    # Warm start the bounding MIP
    convhull_no_good_cuts::Bool                                 # Add no-good cuts to MIP based on the pool solutions

    # Presolving parameters
    presolve_track_time::Bool                                   # Account presolve time for total time usage
    presolve_bt::Bool                                           # Perform bound tightening procedure before the main algorithm (default: true)
    presolve_time_limit::Float64                                # Time limit for presolving (seconds)
    presolve_bt_max_iter::Int                                   # Maximum iterations allowed to perform presolve
    presolve_bt_width_tol::Float64                              # Width tolerance for bound-tightening
    presolve_bt_output_tol::Float64                             # Variable bounds truncation tol (change to precision)
    presolve_bt_algo::Any                                       # Method used for bound tightening procedures, can either be an index of default methods or functional inputs
    presolve_bt_relax_integrality::Bool                                     # Relax the MIP solved in built-in relaxation scheme for time performance
    presolve_bt_mip_time_limit::Float64                         # Time limit for a single MIP solved in the built-in bound tightening algorithm (with partitions)

    # Domain Reduction
    presolve_bp::Bool                                           # Conduct basic bound propagation
    user_parameters::Dict                                       # [INACTIVE] Additional parameters used for user-defined functional inputs

    # Features for Integer Problems (NOTE: no support for int-lin problems)
    int_enable::Bool                                            # Convert integer problem into binary problem
    int_cumulative_disc::Bool                                   #  Cumulatively involve integer variables for discretization

end

function default_options()
        log_level = 1
        time_limit = 1E6
        max_iter = 99
        rel_gap = 1e-4
        gap_ref = :ub
        abs_gap = 1e-6
        tol = 1e-6
        large_bound = 1E6

        nlp_solver = nothing
        minlp_solver = nothing
        mip_solver = nothing

        recognize_convex = true
        bilinear_mccormick = false
        bilinear_convexhull = true
        monomial_convexhull = true

        method_convexification = Array{Function}(undef, 0)
        # method_partition_injection = Array{Function}(undef, 0)
        # term_patterns = Array{Function}(undef, 0)
        # constr_patterns = Array{Function}(undef, 0)

        disc_var_pick = 2                      # By default use the 15-variable selective rule
        disc_ratio = 4
        disc_uniform_rate = 2
        disc_add_partition_method = "adaptive"
        disc_divert_chunks = 5
        disc_abs_width_tol = 1e-4
        disc_rel_width_tol = 1e-6
        disc_consecutive_forbid = false
        disc_ratio_branch = false

        convhull_formulation = "sos2"
        convhull_ebd = false
        convhull_ebd_encode = "default"
        convhull_ebd_ibs = false
        convhull_ebd_link = false
        convhull_warmstart = true
        convhull_no_good_cuts = true

        presolve_track_time = true
        presolve_bt = true
        presolve_time_limit = 900
        presolve_bt_max_iter = 10
        presolve_bt_width_tol = 1e-3
        presolve_bt_output_tol = 1e-5
        presolve_bt_algo = 1
        presolve_bt_relax_integrality = false
        presolve_bt_mip_time_limit = Inf
        presolve_bp = false

        user_parameters = Dict()
        int_enable = false
        int_cumulative_disc = true

    return OptimizerOptions(log_level, time_limit, max_iter, rel_gap, gap_ref, abs_gap, tol, large_bound,
                             nlp_solver, minlp_solver, mip_solver,
                             recognize_convex, bilinear_mccormick, bilinear_convexhull, monomial_convexhull,
                             method_convexification, disc_var_pick, disc_ratio, disc_uniform_rate, disc_add_partition_method, disc_divert_chunks,
                             disc_abs_width_tol, disc_rel_width_tol, disc_consecutive_forbid, disc_ratio_branch,
                             convhull_formulation, convhull_ebd, convhull_ebd_encode, convhull_ebd_ibs, convhull_ebd_link, convhull_warmstart, convhull_no_good_cuts,
                             presolve_track_time, presolve_bt, presolve_time_limit, presolve_bt_max_iter, presolve_bt_width_tol, presolve_bt_output_tol,
                             presolve_bt_algo, presolve_bt_relax_integrality, presolve_bt_mip_time_limit, presolve_bp,
                             user_parameters, int_enable, int_cumulative_disc)
end


mutable struct Optimizer <: MOI.AbstractOptimizer

    options::OptimizerOptions                                   # Options set by user

    # Sub-solver identifier for customized solver option
    nlp_solver_id::AbstractString                               # NLP Solver identifier string
    minlp_solver_id::AbstractString                             # MINLP local solver identifier string
    mip_solver_id::AbstractString                               # MIP solver identifier string

    # User inputs
    num_var_orig::Int                                           # Initial number of variables
    num_cont_var_orig::Int                                      # Initial number of continuous variables
    num_int_var_orig::Int                                       # Initial number of binary/integer variables
    num_constr_orig::Int                                        # Initial number of constraints
    num_lconstr_orig::Int                                       # Initial number of linear constraints
    num_nlconstr_orig::Int                                      # Initial number of non-linear constraints
    var_type_orig::Vector{Symbol}                               # Variable type vector on original variables (only :Bin, :Cont, :Int)
    var_start_orig::Vector{Float64}                             # Variable warm start vector on original variables
    constr_type_orig::Vector{Symbol}                            # Constraint type vector on original variables (only :(==), :(>=), :(<=))
    lin_quad_constraints::Vector{Any}                           # Constraint `func`-in-`set` values
    constr_expr_orig::Vector{Expr}                              # Constraint expressions
    obj_expr_orig::Union{Expr,Number}                           # Objective expression

    # Additional user inputs concerning local solves
    l_var_orig::Vector{Float64}                                 # Variable lower bounds
    u_var_orig::Vector{Float64}                                 # Variable upper bounds
    constraint_bounds_orig::Vector{MOI.NLPBoundsPair}           # Constraint lower bounds
    nl_constraint_bounds_orig::Vector{MOI.NLPBoundsPair}        # Constraint lower bounds
    sense_orig::MOI.OptimizationSense                           # Problem type (:Min, :Max)
    d_orig::Union{Nothing, JuMP.NLPEvaluator}                   # Instance of AbstractNLPEvaluator for evaluating gradient, Hessian-vector products, and Hessians of the Lagrangian
    has_nl_objective::Bool
    objective_function::Union{Nothing, MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}}

    # Additional initial data
    is_obj_linear_orig::Bool                                      # Boolean parameter for type of objective

    # (un-populated options for later use)
    # A_orig::Any                                                 # Linear constraint matrix
    # A_l_orig::Vector{Float64}                                   # Linear constraint matrix LHS
    # A_u_orig::Vector{Float64}                                   # Linear constraint matrix RHS
    # c_orig::Vector{Float64}                                     # Coefficient vector for linear objective
    # num_lconstr_updated::Int                                    # Updated number of linear constraints - includes linear constraints added via @NLconstraint macro
    # num_nlconstr_updated::Int                                   # Updated number of non-linear constraints
    # indices_lconstr_updated::Vector{Int}                        # Indexes of updated linear constraints

    # Local solution model (extra data for each iteration)
    l_var::Vector{Float64}                                      # Updated variable lower bounds for local solve
    u_var::Vector{Float64}                                      # Updated variable upper bounds for local solve

    # MIP bounding model
    model_mip::JuMP.Model                                       # JuMP's convex-MIP model for bounding

    num_var_linear_mip::Int                                     # Number of linear lifted variables required
    num_var_nonlinear_mip::Int                                  # Number of lifted variables
    num_var_disc_mip::Int                                       # Number of variables which are discretized
    num_constr_convex::Int                                      # Number of convex constraints

    # Expression-related and structural property placeholder
    linear_terms::Dict{Any, Any}                                # Dictionary containing details of lifted linear terms
    nonconvex_terms::Dict{Any,Any}                              # Dictionary containing details of lifted non-linear terms
    term_seq::Dict{Int, Any}                                    # Vector-Dictionary for NL terms detection
    nonlinear_constrs::Dict{Any,Any}                            # Dictionary containing details of special constraints
    obj_structure::Symbol                                       # A symbolic indicator of the expression type of objective function
    constr_structure::Vector{Symbol}                            # A vector indicating whether a constraint is with the special structure
    bounding_obj_expr_mip::Union{Expr,Number}                   # Lifted objective expression; if linear, same as obj_expr_orig
    bounding_constr_expr_mip::Vector{Expr}                      # Lifted constraints; if linear, same as corresponding constr_expr_orig
    bounding_obj_mip::Dict{Any, Any}                            # Lifted objective expression in affine form
    bounding_constr_mip::Vector{Dict{Any, Any}}                 # Lifted constraint expressions in affine form

    # Discretization Related options
    candidate_disc_vars::Vector{Int}                            # A vector of all original variable indices that is involved in the nonlinear terms
    discretization::Dict{Any,Any}                               # Discretization points with variable keys
    disc_vars::Vector{Int}                                      # Variables chosen for discretization
    int_vars::Vector{Int}                                       # Index vector of integer variables
    bin_vars::Vector{Int}                                       # Index vector of binary variables

    # Reformulated problem options
    l_var_tight::Vector{Float64}                                # Tightened variable upper bounds
    u_var_tight::Vector{Float64}                                # Tightened variable lower bounds
    var_type::Vector{Symbol}                                    # Updated variable type for local solve

    # Solution information
    presolve_infeasible::Bool                                   # Presolve infeasibility detection flag
    best_bound::Float64                                         # Best bound from MIP
    best_obj::Float64                                           # Best feasible objective value
    initial_warmval::Vector{Float64}                            # Warmstart values set to Alpine
    best_sol::Vector{Float64}                                   # Best feasible solution
    best_bound_sol::Vector{Float64}                             # Best bound solution (arg-min)
    best_rel_gap::Float64                                       # Relative optimality gap = |best_obj - best_bound|/|best_obj|*100
    best_abs_gap::Float64                                       # Absolute gap = |best_obj - best_bound|
    bound_sol_history::Vector{Vector{Float64}}                  # History of bounding solutions limited by "parameter disc_consecutive_forbid"
    bound_sol_pool::Dict{Any, Any}                              # A pool of solutions from solving model_mip

    # Logging information and status
    logs::Dict{Symbol,Any}                                      # Logging information
    detected_feasible_solution::Bool
    detected_bound::Bool
    status::Dict{Symbol, MOI.TerminationStatusCode}             # Detailed status of every iteration in the algorithm
    alpine_status::MOI.TerminationStatusCode                    # Current Alpine's status

    # Constructor for Alpine.Optimizer
    function Optimizer()

        m = new()
        m.options = default_options()
        MOI.empty!(m)

        return m
    end
end

struct NumberOfIterations <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::NumberOfIterations) = true
MOI.get(m::Optimizer, ::NumberOfIterations) = m.logs[:n_iter]

struct NumberOfPresolveIterations <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::NumberOfPresolveIterations) = true
MOI.get(m::Optimizer, ::NumberOfPresolveIterations) = m.logs[:bt_iter]

MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m.alpine_status

function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif m.alpine_status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif m.alpine_status == MOI.LOCALLY_SOLVED
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

MOI.get(::Optimizer, ::MOI.PrimalStatus) = MOI.NO_SOLUTION

MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m.best_obj
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m.best_bound
MOI.get(m::Optimizer, ::MOI.SolveTimeSec) = m.logs[:total_time]

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
    m.nlp_solver_id = ""
    m.minlp_solver_id = ""
    m.mip_solver_id = ""

    m.num_var_orig = 0
    m.num_cont_var_orig = 0
    m.num_int_var_orig = 0
    m.num_constr_orig = 0
    m.num_lconstr_orig = 0
    m.num_nlconstr_orig = 0
    m.var_type_orig = Symbol[]
    m.var_start_orig = Float64[]
    m.constr_type_orig = Symbol[]
    m.lin_quad_constraints = Any[]
    m.constr_expr_orig = Expr[]
    # m.num_lconstr_updated = 0
    # m.num_nlconstr_updated = 0
    # m.indices_lconstr_updated = Int[]

    m.l_var_orig = Float64[]
    m.u_var_orig = Float64[]
    m.constraint_bounds_orig = MOI.NLPBoundsPair[]
    m.nl_constraint_bounds_orig = MOI.NLPBoundsPair[]
    m.sense_orig = MOI.FEASIBILITY_SENSE

    m.d_orig = nothing
    m.has_nl_objective = false
    m.objective_function = nothing

    m.linear_terms = Dict()
    m.nonconvex_terms = Dict()
    m.term_seq = Dict()
    m.nonlinear_constrs = Dict()
    m.candidate_disc_vars = Int[]
    m.bounding_constr_expr_mip = []
    m.bounding_constr_mip = []
    m.disc_vars = []
    m.int_vars = []
    m.bin_vars = []
    m.discretization = Dict()
    m.num_var_linear_mip = 0
    m.num_var_nonlinear_mip = 0
    m.num_var_disc_mip = 0
    m.num_constr_convex = 0
    m.constr_structure = Symbol[]
    m.best_bound_sol = []
    m.bound_sol_history = []
    m.presolve_infeasible = false
    m.bound_sol_history = Vector{Vector{Float64}}(undef, m.options.disc_consecutive_forbid)

    m.best_obj = Inf
    m.initial_warmval = Float64[]
    m.best_sol = Float64[]
    m.best_bound = -Inf
    m.best_rel_gap = Inf
    m.best_abs_gap = Inf
    m.alpine_status = MOI.OPTIMIZE_NOT_CALLED

    create_status!(m)
    create_logs!(m)
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOIU.default_copy_to(model, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Alpine"

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

_lower(set::MOI.EqualTo) = set.value
_upper(set::MOI.EqualTo) = set.value
_lower(set::MOI.LessThan) = nothing
_upper(set::MOI.LessThan) = set.upper
_lower(set::MOI.GreaterThan) = set.lower
_upper(set::MOI.GreaterThan) = nothing
_lower(set::MOI.Interval) = set.lower
_upper(set::MOI.Interval) = set.upper

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
    if is_max_sense(model)
        model.best_obj = -Inf
        model.best_bound = Inf
    elseif is_min_sense(model)
        model.best_obj = Inf
        model.best_bound = -Inf
    else
        error("Feasibility sense not supported yet by Alpine.")
    end
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where F
    model.objective_function = func
end

function MOI.set(m::Optimizer, ::MOI.NLPBlock, block)
    m.d_orig = block.evaluator
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

function load!(m::Optimizer)
    # Initialize NLP interface
    MOI.initialize(m.d_orig, [:Grad, :Jac, :Hess, :HessVec, :ExprGraph]) # Safety scheme for sub-solvers re-initializing the NLPEvaluator

    # Collect objective & constraint expressions
    if m.has_nl_objective
        m.obj_expr_orig = expr_isolate_const(_variable_index_to_index(MOI.objective_expr(m.d_orig))) # see in nlexpr.jl if this expr isolation has any issue
    elseif m.objective_function isa Nothing
        m.obj_expr_orig = Expr(:call, :+)
    else
        m.obj_expr_orig = _moi_function_to_expr(m.objective_function)
    end

    # Collect original variable type and build dynamic variable type space
    m.var_type = copy(m.var_type_orig)
    m.int_vars = [i for i in 1:m.num_var_orig if m.var_type[i] == :Int]
    m.bin_vars = [i for i in 1:m.num_var_orig if m.var_type[i] == :Bin]

    if !isempty(m.int_vars) || !isempty(m.bin_vars)
        (Alp.get_option(m, :minlp_solver) === nothing) && (error("No MINLP local solver specified; use option 'minlp_solver' to specify a MINLP local solver"))
    end

    m.num_constr_orig += length(m.nl_constraint_bounds_orig)
    m.num_nlconstr_orig += length(m.nl_constraint_bounds_orig)
    append!(m.constraint_bounds_orig, m.nl_constraint_bounds_orig)
    for i in eachindex(m.nl_constraint_bounds_orig)
        push!(m.constr_expr_orig, _variable_index_to_index(MOI.constraint_expr(m.d_orig, i)))
        push!(m.constr_structure, :generic_nonlinear)
    end

    # Summarize constraints information in original model
    m.constr_type_orig = Array{Symbol}(undef, m.num_constr_orig)

    for i in 1:m.num_constr_orig
        if m.constraint_bounds_orig[i].lower > -Inf && m.constraint_bounds_orig[i].upper < Inf
            m.constr_type_orig[i] = :(==)
        elseif m.constraint_bounds_orig[i].lower > -Inf
            m.constr_type_orig[i] = :(>=)
        else
            m.constr_type_orig[i] = :(<=)
        end
    end

    # Initialize recognizable structure properties with :none
    m.obj_structure = :none

    @assert m.num_constr_orig == m.num_nlconstr_orig + m.num_lconstr_orig
    m.is_obj_linear_orig = !m.has_nl_objective && m.objective_function isa MOI.ScalarAffineFunction{Float64}
    m.is_obj_linear_orig ? (m.obj_structure = :generic_linear) : (m.obj_structure = :generic_nonlinear)
    isa(m.obj_expr_orig, Number) && (m.obj_structure = :constant)

    # populate data to create the bounding model
    recategorize_var(m)             # Initial round of variable re-categorization

    :Int in m.var_type_orig && error("Alpine does not support MINLPs with generic integer (non-binary) variables yet! Try Juniper.jl for finding a local feasible solution")
    :Int in m.var_type_orig ? Alp.set_option(m, :int_enable, true) : Alp.set_option(m, :int_enable, false) # Separator for safer runs

    # Conduct solver-dependent detection
    fetch_mip_solver_identifier(m)
    (Alp.get_option(m, :nlp_solver) !== nothing) && (fetch_nlp_solver_identifier(m))
    (Alp.get_option(m, :minlp_solver) !== nothing) && (fetch_minlp_solver_identifier(m))

    # Solver Dependent Options
    if m.mip_solver_id != :Gurobi
        Alp.get_option(m, :convhull_warmstart) == false
        Alp.get_option(m, :convhull_no_good_cuts) == false
    end

    # Main Algorithmic Initialization
    Alp.process_expr(m)                         # Compact process of every expression
    Alp.init_tight_bound(m)                     # Initialize bounds for algorithmic processes
    Alp.resolve_var_bounds(m)                   # resolve lifted var bounds
    Alp.pick_disc_vars(m)                       # Picking variables to be discretized
    Alp.init_disc(m)                            # Initialize discretization dictionaries

    # Turn-on bt presolve if variables are not discrete
    if isempty(m.int_vars) && length(m.bin_vars) <= 50 && m.num_var_orig <= 10000 && length(m.candidate_disc_vars)<=300 && Alp.get_option(m, :presolve_bt) == nothing
        Alp.set_option(m, :presolve_bt, true)
        println("Automatically turning on bound-tightening presolve")
    elseif Alp.get_option(m, :presolve_bt) == nothing  # If no use indication
        Alp.set_option(m, :presolve_bt, false)
    end

    if length(m.bin_vars) > 200 || m.num_var_orig > 2000
        println("Automatically turning OFF 'disc_ratio_branch' due to the size of the problem")
        Alp.set_option(m, :disc_ratio_branch, false)
    end

    # Initialize the solution pool
    m.bound_sol_pool = initialize_solution_pool(m, 0)  # Initialize the solution pool

    # Check if any illegal term exist in the warm-solution
    any(isnan, m.best_sol) && (m.best_sol = zeros(length(m.best_sol)))

    # Initialize log
    logging_summary(m)

    return
end


function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)

    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    return model.best_sol[vi.value]

end

MOI.get(model::Optimizer, ::MOI.ResultCount) = model.alpine_status == MOI.OPTIMIZE_NOT_CALLED ? 0 : 1

MOI.is_valid(model::Alpine.Optimizer, vi::MOI.VariableIndex) = 1 <= vi.value <= model.num_var_orig
