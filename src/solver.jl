export AlpineSolver

mutable struct Optimizer <: MOI.AbstractOptimizer

    # Parameters for tuning Alpine

    # basic solver parameters
    loglevel::Int                                               # Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info
    timeout::Float64                                            # Time limit for algorithm (in seconds)
    maxiter::Int                                                # Target Maximum Iterations
    relgap::Float64                                             # Relative optimality gap termination condition
    gapref::Symbol                                              # Relative gap reference point (options: [:ub, :lb])
    absgap::Float64                                             # Absolute optimality gap termination condition
    tol::Float64                                                # Numerical tol used in the algorithmic process
    largebound::Float64                                         # Large bounds for problems with unbounded variables

    # convexification method tuning
    recognize_convex::Bool                                      # Recognize convex expressions in parsing objective functions and constraints
    bilinear_mccormick::Bool                                    # Convexify bilinear terms using piecwise McCormick representation
    bilinear_convexhull::Bool                                   # Convexify bilinear terms using lambda representation
    monomial_convexhull::Bool                                   # Convexify monomial terms using convex-hull representation

    # expression-based user-inputs
    method_convexification::Array{Function}                     # Array of functions that user can choose to convexify specific non-linear terms : no over-ride privilege
    method_partition_injection::Array{Function}                 # Array of functions for special methods to add partitions to variables under complex conditions
    term_patterns::Array{Function}                              # Array of functions that user can choose to parse/recognize nonlinear terms in constraint expressions
    constr_patterns::Array{Function}                            # Array of functions that user can choose to parse/recognize structural constraint from expressions

    # parameters used in the partitioning algorithm
    disc_ratio::Any                                             # Discretization ratio parameter (use a fixed value for now, later switch to a function)
    disc_uniform_rate::Int                                      # Discretization rate parameter when using uniform partitions
    disc_var_pick::Any                                          # Algorithm for choosing the variables to discretize: 1 for minimum vertex cover, 0 for all variables
    disc_divert_chunks::Int                                     # How many uniform partitions to construct
    disc_add_partition_method::Any                              # Additional methods to add discretization
    disc_abs_width_tol::Float64                                 # Absolute tolerance used when setting up partition/discretization
    disc_rel_width_tol::Float64                                 # Relative width tolerance when setting up partition/discretization
    disc_consecutive_forbid::Int                                # Prevent bounding model to add partitions consecutively in the same region when bounds do not improve
    disc_ratio_branch::Bool                                     # Branching tests for picking fixed the discretization ratio

    # MIP Formulation Parameters
    convhull_formulation::String                                # MIP Formulation for the relaxation
    convhull_warmstart::Bool                                    # Warm start the bounding MIP
    convhull_no_good_cuts::Bool                                 # Add no-good cuts to MIP based on the pool solutions
    convhull_ebd::Bool                                          # Enable embedding formulation
    convhull_ebd_encode::Any                                    # Encoding method used for convhull_ebd
    convhull_ebd_ibs::Bool                                      # Enable independent branching scheme
    convhull_ebd_link::Bool                                     # Linking constraints between x and α, type 1 uses hierarchical and type 2 uses big-m

    # Presolving Parameters
    presolve_track_time::Bool                                   # Account presolve time for total time usage
    presolve_bt::Bool                                           # Perform bound tightening procedure before the main algorithm (default: true)
    presolve_timeout::Float64                                   # Time limit for presolving (seconds)
    presolve_maxiter::Int                                       # Maximum iterations allowed to perform presolve (vague in parallel mode)
    presolve_bt_width_tol::Float64                              # Width tolerance for bound-tightening
    presolve_bt_output_tol::Float64                             # Variable bounds truncation tol (change to precision)
    presolve_bt_algo::Any                                       # Method used for bound tightening procedures, can either be an index of default methods or functional inputs
    presolve_bt_relax::Bool                                     # Relax the MIP solved in built-in relaxation scheme for time performance
    presolve_bt_mip_timeout::Float64                            # Regulate the time limit for a single MIP solved in the built-in bound tightening algorithm

    # Domain Reduction
    presolve_bp::Bool                                           # Conduct basic bound propagation
    presolve_infeasible::Bool                                   # Presolve infeasibility detection flag
    user_parameters::Dict                                       # Additional parameters used for user-defined functional inputs

    # Features for Integer Problems (NOTE: no support for int-lin problems)
    int_enable::Bool                                            # Convert integer problem into binary problem
    int_cumulative_disc::Bool                                   # [INACTIVE] Cumulatively involve integer variables for discretization
    int_fully_disc::Bool                                        # [INACTIVE] Construct equivalent formulation for integer variables

    # add all the solver options
    nlp_solver                                                  # Local continuous NLP solver for solving NLPs at each iteration
    minlp_solver                                                # Local MINLP solver for solving MINLPs at each iteration
    mip_solver                                                  # MIP solver for successive lower bound solves

    # Sub-solver identifier for customized solver option
    nlp_solver_id::AbstractString                               # NLP Solver identifier string
    minlp_solver_id::AbstractString                             # MINLP local solver identifier string
    mip_solver_id::AbstractString                               # MIP solver identifier string

    # user provided inputs
    num_var_orig::Int                                           # Initial number of variables
    num_cont_var_orig::Int                                      # Initial number of continuous variables
    num_int_var_orig::Int                                       # Initial number of binary/integer variables
    num_constr_orig::Int                                        # Initial number of constraints
    num_lconstr_orig::Int                                       # Initial number of linear constraints
    num_nlconstr_orig::Int                                      # Initial number of non-linear constraints
    var_type_orig::Vector{Symbol}                               # Variable type vector on original variables (only :Bin, :Cont, :Int)
    var_start_orig::Vector{Float64}                             # Variable warm start vector on original variables
    constr_type_orig::Vector{Symbol}                            # Constraint type vector on original variables (only :(==), :(>=), :(<=))
    constr_expr_orig::Vector{Expr}                              # Constraint expressions
    obj_expr_orig::Union{Expr,Number}                                         # Objective expression

    # additional user inputs useful for local solves
    l_var_orig::Vector{Float64}                                 # Variable lower bounds
    u_var_orig::Vector{Float64}                                 # Variable upper bounds
    l_constr_orig::Vector{Float64}                              # Constraint lower bounds
    u_constr_orig::Vector{Float64}                              # Constraint upper bounds
    sense_orig::Symbol                                          # Problem type (:Min, :Max)
    d_orig::JuMP.NLPEvaluator                                   # Instance of AbstractNLPEvaluator for evaluating gradient, Hessian-vector products, and Hessians of the Lagrangian

    # additional initial data that may be useful later (not populated)
    A_orig::Any                                                 # Linear constraint matrix
    A_l_orig::Vector{Float64}                                   # Linear constraint matrix LHS
    A_u_orig::Vector{Float64}                                   # Linear constraint matrix RHS
    is_obj_linear_orig::Bool                                    # Bool variable for type of objective
    c_orig::Vector{Float64}                                     # Coefficient vector for linear objective
    num_lconstr_updated::Int                                    # Updated number of linear constraints - includes linear constraints added via @NLconstraint macro
    num_nlconstr_updated::Int                                   # Updated number of non-linear constraints
    indexes_lconstr_updated::Vector{Int}                        # Indexes of updated linear constraints

    # local solution model extra data for each iteration
    l_var::Vector{Float64}                                      # Updated variable lower bounds for local solve
    u_var::Vector{Float64}                                      # Updated variable upper bounds for local solve

    # MIP bounding model
    model_mip::JuMP.Model                                       # JuMP's convex-MIP model for bounding

    num_var_linear_mip::Int                                     # Number of linear lifted variables required
    num_var_nonlinear_mip::Int                                  # Number of lifted variables
    num_var_disc_mip::Int                                       # Number of variables which are discretized
    num_constr_convex::Int                                      # Number of convex constraints

    # expression-related and structural property placeholder
    linear_terms::Dict{Any, Any}                                # Dictionary containing details of lifted linear terms
    nonconvex_terms::Dict{Any,Any}                              # Dictionary containing details of lifted non-linear terms
    term_seq::Dict{Int, Any}                                    # Vector-Dictionary for NL terms detection
    nonlinear_constrs::Dict{Any,Any}                            # Dictionary containing details of special constraints
    obj_structure::Symbol                                       # A symbolic indicator of the expression type of objective function
    constr_structure::Vector{Symbol}                            # A vector indicating whether a constraint is with the special structure
    bounding_obj_expr_mip::Union{Expr,Number}                                 # Lifted objective expression; if linear, same as obj_expr_orig
    bounding_constr_expr_mip::Vector{Expr}                      # Lifted constraints; if linear, same as corresponding constr_expr_orig
    bounding_obj_mip::Dict{Any, Any}                            # Lifted objective expression in affine form
    bounding_constr_mip::Vector{Dict{Any, Any}}                 # Lifted constraint expressions in affine form

    # Discretization Related
    candidate_disc_vars::Vector{Int}                            # A vector of all original variable indices that is involved in the nonlinear terms
    discretization::Dict{Any,Any}                               # Discretization points with variable keys
    disc_vars::Vector{Int}                                      # Variables chosen for discretization
    int_vars::Vector{Int}                                       # Index vector of integer variables
    bin_vars::Vector{Int}                                       # Index vector of binary variables

    # Reformulated problem
    l_var_tight::Vector{Float64}                                # Tightened variable upper bounds
    u_var_tight::Vector{Float64}                                # Tightened variable lower bounds
    var_type::Vector{Symbol}                                    # Updated variable type for local solve

    # Solution information
    best_bound::Float64                                         # Best bound from MIP
    best_obj::Float64                                           # Best feasible objective value
    best_sol::Vector{Float64}                                   # Best feasible solution
    best_bound_sol::Vector{Float64}                             # Best bound solution (arg-min)
    best_rel_gap::Float64                                       # Relative optimality gap = |best_obj - best_bound|/|best_obj|*100
    best_abs_gap::Float64                                       # Absolute gap = |best_obj - best_bound|
    bound_sol_history::Vector{Vector{Float64}}                  # History of bounding solutions limited by "parameter disc_consecutive_forbid"
    bound_sol_pool::Dict{Any, Any}                              # A pool of solutions from solving model_mip

    # Logging information and status
    logs::Dict{Symbol,Any}                                      # Logging information
    status::Dict{Symbol,Symbol}                                 # Detailed status of every iteration in the algorithm
    alpine_status::Symbol                                       # Current Alpine's status

    # constructor
    function Optimizer()

        m = new()

        m.loglevel = 1
        m.timeout = Inf
        m.maxiter = 99
        m.relgap = 1e-4
        m.gapref = :ub
        m.absgap = 1e-6
        m.tol = 1e-6
        m.largebound = 1e4

        m.nlp_solver = nothing
        m.minlp_solver = nothing
        m.mip_solver = nothing

        m.recognize_convex = true
        m.bilinear_mccormick = false
        m.bilinear_convexhull = true
        m.monomial_convexhull = true

        m.method_convexification = Array{Function}(undef, 0)
        m.method_partition_injection = Array{Function}(undef, 0)
        m.term_patterns = Array{Function}(undef, 0)
        m.constr_patterns = Array{Function}(undef, 0)

        m.disc_var_pick = 2                      # By default use the 15-variable selective rule
        m.disc_ratio = 4
        m.disc_uniform_rate = 2
        m.disc_add_partition_method = "adaptive"
        m.disc_divert_chunks = 5
        m.disc_abs_width_tol = 1e-4
        m.disc_rel_width_tol = 1e-6
        m.disc_consecutive_forbid = 0
        m.disc_ratio_branch=false

        m.convhull_formulation = "sos2"
        m.convhull_ebd = false
        m.convhull_ebd_encode = "default"
        m.convhull_ebd_ibs = false
        m.convhull_ebd_link = false
        m.convhull_warmstart = true
        m.convhull_no_good_cuts = true

        m.presolve_track_time = true
        m.presolve_bt = true
        m.presolve_timeout = 900
        m.presolve_maxiter = 10
        m.presolve_bt_width_tol = 1e-3
        m.presolve_bt_output_tol = 1e-5
        m.presolve_bt_algo = 1
        m.presolve_bt_relax = false
        m.presolve_bt_mip_timeout = Inf

        m.presolve_bp = true

        m.user_parameters = Dict()
        m.int_enable = false
        m.int_cumulative_disc = true
        m.int_fully_disc = false

        m.num_var_orig = 0
        m.num_cont_var_orig = 0
        m.num_int_var_orig = 0
        m.num_constr_orig = 0
        m.num_lconstr_orig = 0
        m.num_nlconstr_orig = 0
        m.var_type_orig = Symbol[]
        m.var_start_orig = Float64[]
        m.constr_type_orig = Symbol[]
        m.constr_expr_orig = Expr[]
        m.num_lconstr_updated = 0
        m.num_nlconstr_updated = 0
        m.indexes_lconstr_updated = Int[]

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
        m.constr_structure = []
        m.best_bound_sol = []
        m.bound_sol_history = []
        m.presolve_infeasible = false
        m.bound_sol_history = Vector{Vector{Float64}}(undef, m.disc_consecutive_forbid)

        m.best_obj = Inf
        m.best_bound = -Inf
        m.best_rel_gap = Inf
        m.best_abs_gap = Inf
        m.alpine_status = :NotLoaded

        create_status!(m)
        create_logs!(m)

        return m
    end
end

function MOI.is_empty(model::Optimizer)
    return iszero(model.num_var_orig)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    setproperty!(model, Symbol(param.name), value)
end

function MOI.add_variables(model::Optimizer, n::Int)
    return [MOI.add_variable(model) for i in 1:n]
end
function MOI.add_variable(model::Optimizer)
    model.num_var_orig += 1
    push!(model.l_var_orig, -Inf)
    push!(model.u_var_orig, Inf)
    return MOI.VariableIndex(model.num_var_orig)
end

const SCALAR_SET = Union{MOI.EqualTo{Float64}, MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.Interval{Float64}}

_lower(set::MOI.EqualTo) = set.value
_upper(set::MOI.EqualTo) = set.value
_lower(set::MOI.LessThan) = nothing
_upper(set::MOI.LessThan) = set.upper
_lower(set::MOI.GreaterThan) = set.lower
_upper(set::MOI.GreaterThan) = nothing
_lower(set::MOI.Interval) = set.lower
_upper(set::MOI.Interval) = set.upper

function MOI.add_constraint(model::Optimizer, f::MOI.SingleVariable, set::SCALAR_SET)
    vi = f.variable
    l = _lower(set)
    if l !== nothing
        model.l_var_orig[vi.value] = l
    end
    u = _upper(set)
    if u !== nothing
        model.u_var_orig[vi.value] = u
    end
    return MOI.ConstraintIndex{typeof(f), typeof(set)}(vi.value)
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense)
    if sense == MOI.MAX_SENSE
        model.sense_orig = :Max
        m.best_obj = -Inf
        m.best_bound = Inf
    else
        model.sense_orig = :Min
        m.best_obj = Inf
        m.best_bound = -Inf
    end
end
function MOI.set(model::Optimizer, ::MOI.NLPBlock, block)
    m.d_orig = block.evaluator
    m.num_constr_orig = length(block.constraint_bounds)
    m.l_constr_orig = [p.lower for p in block.constraint_bounds]
    m.u_constr_orig = [p.upper for p in block.constraint_bounds]
    return
end

function MOI.optimize!(m::Optimizer)
    # Initialize NLP interface
    interface_init_nonlinear_data(m.d_orig)

    # Collect objective & constraint expressions
    m.obj_expr_orig = expr_isolate_const(interface_get_obj_expr(m.d_orig)) # see in nlexpr.jl if this expr isolation has any issue

    for i in 1:m.num_constr_orig
        push!(m.constr_expr_orig, interface_get_constr_expr(m.d_orig, i))
    end

    # Collect original variable type and build dynamic variable type space
    m.var_type_orig = [getcategory(Variable(d.m, i)) for i in 1:m.num_var_orig]
    m.var_type = copy(m.var_type_orig)
    m.int_vars = [i for i in 1:m.num_var_orig if m.var_type[i] == :Int]
    m.bin_vars = [i for i in 1:m.num_var_orig if m.var_type[i] == :Bin]

    if !isempty(m.int_vars) || !isempty(m.bin_vars)
        (m.minlp_solver === nothing) && (error("No MINLP local solver specified; use minlp_solver to specify a MINLP local solver"))
    end

    # Summarize constraints information in original model
    m.constr_type_orig = Array{Symbol}(undef, m.num_constr_orig)

    for i in 1:m.num_constr_orig
        if l_constr[i] > -Inf && u_constr[i] < Inf
            m.constr_type_orig[i] = :(==)
        elseif l_constr[i] > -Inf
            m.constr_type_orig[i] = :(>=)
        else
            m.constr_type_orig[i] = :(<=)
        end
    end

    # Initialize recognizable structure properties with :none
    m.obj_structure = :none
    m.constr_structure = [:none for i in 1:m.num_constr_orig]
    for i = 1:m.num_constr_orig
        if interface_is_constr_linear(m.d_orig, i)
            m.num_lconstr_orig += 1
            m.constr_structure[i] = :generic_linear
        else
            m.num_nlconstr_orig += 1
            m.constr_structure[i] = :generic_nonlinear
        end
    end

    @assert m.num_constr_orig == m.num_nlconstr_orig + m.num_lconstr_orig
    m.is_obj_linear_orig = interface_is_obj_linear(m.d_orig)
    m.is_obj_linear_orig ? (m.obj_structure = :generic_linear) : (m.obj_structure = :generic_nonlinear)
    isa(m.obj_expr_orig, Number) && (m.obj_structure = :constant)

    # populate data to create the bounding model
    recategorize_var(m)             # Initial round of variable re-categorization

    :Int in m.var_type_orig && @warn "Alpine's support for integer variables is experimental"
    :Int in m.var_type_orig ? m.int_enable = true : m.int_enable = false # Separator for safer runs

    # Conduct solver-dependent detection
    fetch_mip_solver_identifier(m)
    (m.nlp_solver != empty_solver) && (fetch_nlp_solver_identifier(m))
    (m.minlp_solver != empty_solver) && (fetch_minlp_solver_identifier(m))

    # Solver Dependent Options
    if m.mip_solver_id != :Gurobi
        m.convhull_warmstart == false
        m.convhull_no_good_cuts == false
    end

    # Main Algorithmic Initialization
    process_expr(m)                         # Compact process of every expression
    init_tight_bound(m)                     # Initialize bounds for algorithmic processes
    resolve_var_bounds(m)                   # resolve lifted var bounds
    pick_disc_vars(m)                       # Picking variables to be discretized
    init_disc(m)                            # Initialize discretization dictionaries

    # Turn-on bt presolver if variables are not discrete
    if isempty(m.int_vars) && length(m.bin_vars) <= 50 && m.num_var_orig <= 10000 && length(m.candidate_disc_vars)<=300 && m.presolve_bt == nothing
        m.presolve_bt = true
        println("Automatically turning on bound-tightening presolver...")
    elseif m.presolve_bt == nothing  # If no use indication
        m.presolve_bt = false
    end

    if length(m.bin_vars) > 200 || m.num_var_orig > 2000
        println("Automatically turning OFF ratio branching due to the size of the problem")
        m.disc_ratio_branch=false
    end

    # Initialize the solution pool
    m.bound_sol_pool = initialize_solution_pool(m, 0)  # Initialize the solution pool

    # Record the initial solution from the warm-starting value, if any
    m.best_sol = m.d_orig.m.colVal

    # Check if any illegal term exist in the warm-solution
    any(isnan, m.best_sol) && (m.best_sol = zeros(length(m.best_sol)))

    # Initialize log
    logging_summary(m)

    optimize!(m)

    return
end
