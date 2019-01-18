export PODSolver

mutable struct PODNonlinearModel <: MathProgBase.AbstractNonlinearModel

    # external developer parameters for testing and debugging
    colorful_pod::Any                                           # Colorful output (remove)

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
    nlp_solver::MathProgBase.AbstractMathProgSolver             # Local continuous NLP solver for solving NLPs at each iteration
    minlp_solver::MathProgBase.AbstractMathProgSolver           # Local MINLP solver for solving MINLPs at each iteration
    mip_solver::MathProgBase.AbstractMathProgSolver             # MIP solver for successive lower bound solves

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
    obj_expr_orig::Expr                                         # Objective expression

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
    bounding_obj_expr_mip::Expr                                 # Lifted objective expression; if linear, same as obj_expr_orig
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
    pod_status::Symbol                                          # Current Alpine's status

    # constructor
    function PODNonlinearModel(colorful_pod,
                                loglevel, timeout, maxiter, relgap, gapref, absgap, tol, largebound,
                                nlp_solver,
                                minlp_solver,
                                mip_solver,
                                recognize_convex,
                                bilinear_mccormick,
                                bilinear_convexhull,
                                monomial_convexhull,
                                method_convexification,
                                method_partition_injection,
                                term_patterns,
                                constr_patterns,
                                disc_var_pick,
                                disc_ratio,
                                disc_uniform_rate,
                                disc_add_partition_method,
                                disc_divert_chunks,
                                disc_abs_width_tol,
                                disc_rel_width_tol,
                                disc_consecutive_forbid,
                                disc_ratio_branch,
                                convhull_formulation,
                                convhull_ebd,
                                convhull_ebd_encode,
                                convhull_ebd_ibs,
                                convhull_ebd_link,
                                convhull_warmstart,
                                convhull_no_good_cuts,
                                presolve_track_time,
                                presolve_bt,
                                presolve_timeout,
                                presolve_maxiter,
                                presolve_bt_width_tol,
                                presolve_bt_output_tol,
                                presolve_bt_algo,
                                presolve_bt_relax,
                                presolve_bt_mip_timeout,
                                presolve_bp,
                                user_parameters,
                                int_enable,
                                int_cumulative_disc,
                                int_fully_disc)

        m = new()

        m.colorful_pod = colorful_pod

        m.loglevel = loglevel
        m.timeout = timeout
        m.maxiter = maxiter
        m.relgap = relgap
        m.gapref = gapref
        m.absgap = absgap
        m.tol = tol
        m.largebound = largebound

        m.recognize_convex = recognize_convex
        m.bilinear_mccormick = bilinear_mccormick
        m.bilinear_convexhull = bilinear_convexhull
        m.monomial_convexhull = monomial_convexhull

        m.method_convexification = method_convexification
        m.method_partition_injection = method_partition_injection
        m.term_patterns = term_patterns
        m.constr_patterns = constr_patterns

        m.disc_var_pick = disc_var_pick
        m.disc_ratio = disc_ratio
        m.disc_uniform_rate = disc_uniform_rate
        m.disc_add_partition_method = disc_add_partition_method
        m.disc_divert_chunks = disc_divert_chunks
        m.disc_abs_width_tol = disc_abs_width_tol
        m.disc_rel_width_tol = disc_rel_width_tol
        m.disc_consecutive_forbid = disc_consecutive_forbid
        m.disc_ratio_branch = disc_ratio_branch

        m.convhull_formulation = convhull_formulation
        m.convhull_ebd = convhull_ebd
        m.convhull_ebd_encode = convhull_ebd_encode
        m.convhull_ebd_ibs = convhull_ebd_ibs
        m.convhull_ebd_link = convhull_ebd_link
        m.convhull_warmstart = convhull_warmstart
        m.convhull_no_good_cuts = convhull_no_good_cuts

        m.presolve_track_time = presolve_track_time
        m.presolve_bt = presolve_bt
        m.presolve_timeout = presolve_timeout
        m.presolve_maxiter = presolve_maxiter
        m.presolve_bt_width_tol = presolve_bt_width_tol
        m.presolve_bt_output_tol = presolve_bt_output_tol
        m.presolve_bt_algo = presolve_bt_algo
        m.presolve_bt_relax = presolve_bt_relax
        m.presolve_bt_mip_timeout = presolve_bt_mip_timeout

        m.presolve_bp = presolve_bp

        m.nlp_solver = nlp_solver
        m.minlp_solver = minlp_solver
        m.mip_solver = mip_solver

        m.user_parameters = user_parameters
        m.int_enable = int_enable
        m.int_cumulative_disc = int_cumulative_disc
        m.int_fully_disc = int_fully_disc

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
        m.pod_status = :NotLoaded

        create_status!(m)
        create_logs!(m)

        return m
    end
end

mutable struct UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

const empty_solver = UnsetSolver()

mutable struct PODSolver <: MathProgBase.AbstractMathProgSolver

    colorful_pod::Any

    loglevel::Int
    timeout::Float64
    maxiter::Int
    relgap::Float64
    gapref::Symbol
    absgap::Float64
    tol::Float64
    largebound::Float64

    nlp_solver::MathProgBase.AbstractMathProgSolver
    minlp_solver::MathProgBase.AbstractMathProgSolver
    mip_solver::MathProgBase.AbstractMathProgSolver

    recognize_convex::Bool
    bilinear_mccormick::Bool
    bilinear_convexhull::Bool
    monomial_convexhull::Bool

    method_convexification::Array{Function}
    method_partition_injection::Array{Function}
    term_patterns::Array{Function}
    constr_patterns::Array{Function}

    disc_var_pick::Any
    disc_ratio::Any
    disc_uniform_rate::Int
    disc_add_partition_method::Any
    disc_divert_chunks::Int
    disc_abs_width_tol::Float64
    disc_rel_width_tol::Float64
    disc_consecutive_forbid::Int
    disc_ratio_branch::Bool

    convhull_formulation::String
    convhull_ebd::Bool
    convhull_ebd_encode::Any
    convhull_ebd_ibs::Bool
    convhull_ebd_link::Bool
    convhull_warmstart::Bool
    convhull_no_good_cuts::Bool

    presolve_track_time::Bool
    presolve_bt::Bool
    presolve_timeout::Float64
    presolve_maxiter::Int
    presolve_bt_width_tol::Float64
    presolve_bt_output_tol::Float64
    presolve_bt_algo::Any
    presolve_bt_relax::Bool
    presolve_bt_mip_timeout::Float64

    presolve_bp::Bool

    user_parameters::Dict
    int_enable::Bool
    int_cumulative_disc::Bool
    int_fully_disc::Bool

    # other options to be added later on
end

function PODSolver(;
    colorful_pod = false,

    loglevel = 1,
    timeout = Inf,
    maxiter = 99,
    relgap = 1e-4,
    gapref = :ub,
    absgap = 1e-6,
    tol = 1e-6,
    largebound = 1e4,

    nlp_solver = empty_solver,
    minlp_solver = empty_solver,
    mip_solver = empty_solver,

    recognize_convex = true,
    bilinear_mccormick = false,
    bilinear_convexhull = true,
    monomial_convexhull = true,

    method_convexification = Array{Function}(undef, 0),
    method_partition_injection = Array{Function}(undef, 0),
    term_patterns = Array{Function}(undef, 0),
    constr_patterns = Array{Function}(undef, 0),

    disc_var_pick = 2,                      # By default use the 15-variable selective rule
    disc_ratio = 4,
    disc_uniform_rate = 2,
    disc_add_partition_method = "adaptive",
    disc_divert_chunks = 5,
    disc_abs_width_tol = 1e-4,
    disc_rel_width_tol = 1e-6,
    disc_consecutive_forbid = 0,
    disc_ratio_branch=false,

    convhull_formulation = "sos2",
    convhull_ebd = false,
    convhull_ebd_encode = "default",
    convhull_ebd_ibs = false,
    convhull_ebd_link = false,
    convhull_warmstart = true,
    convhull_no_good_cuts = true,

    presolve_track_time = true,
    presolve_maxiter = 10,
    presolve_bt = true,
    presolve_timeout = 900,
    presolve_bt_width_tol = 1e-3,
    presolve_bt_output_tol = 1e-5,
    presolve_bt_algo = 1,
    presolve_bt_relax = false,
    presolve_bt_mip_timeout = Inf,

    presolve_bp = true,

    user_parameters = Dict(),
    int_enable = false,
    int_cumulative_disc = true,
    int_fully_disc = false,

    kwargs...
    )

    # Keyword arguments screening
    unsupported_kwargs = Dict(kwargs)
    !isempty(keys(unsupported_kwargs)) && error("Detected unsupported keyword arguments: $(keys(unsupport_opts))")

    if nlp_solver == empty_solver && minlp_solver == empty_solver
        error("No NLP and MINLP local solver specified; use nlp_solver or minlp_solver to set a local solver dependening on the problem type\n")
    end 
    mip_solver == empty_solver && error("No MIP solver specififed; use mip_solver to set a mip solver\n")

    if nlp_solver != empty_solver && !applicable(MathProgBase.NonlinearModel, nlp_solver)
        error("NLP local solver $(s.nlp_solver) is not supported by JuMP; use a JuMP-supoorted NLP local solver\n")
    end 

    if minlp_solver != empty_solver && !applicable(MathProgBase.NonlinearModel, minlp_solver)
        error("MINLP local solver $(s.minlp_solver) is not supported by JuMP; use a JuMP-supoorted MINLP local solver\n")
    end 

    gapref in [:ub, :lb] || error("the option gapref can take a value only in [:ub, :lb]")

    # String Code Conversion (have to change this to consistently use Symbols and not give user many options)
    if disc_var_pick in ["ncvar_collect_nodes", "all", "max"]
        disc_var_pick = 0
    elseif disc_var_pick in ["min_vertex_cover","min"]
        disc_var_pick = 1
    elseif disc_var_pick == "selective"
        disc_var_pick = 2
    elseif disc_var_pick == "dynamic"
        disc_var_pick = 3
    end

    # Deep-copy the solvers because we may change option values inside Alpine
    PODSolver(colorful_pod,
        loglevel, timeout, maxiter, relgap, gapref, absgap, tol, largebound,
        deepcopy(nlp_solver),
        deepcopy(minlp_solver),
        deepcopy(mip_solver),
        recognize_convex,
        bilinear_mccormick,
        bilinear_convexhull,
        monomial_convexhull,
        method_convexification,
        method_partition_injection,
        term_patterns,
        constr_patterns,
        disc_var_pick,
        disc_ratio,
        disc_uniform_rate,
        disc_add_partition_method,
        disc_divert_chunks,
        disc_abs_width_tol,
        disc_rel_width_tol,
        disc_consecutive_forbid,
        disc_ratio_branch,
        convhull_formulation,
        convhull_ebd,
        convhull_ebd_encode,
        convhull_ebd_ibs,
        convhull_ebd_link,
        convhull_warmstart,
        convhull_no_good_cuts,
        presolve_track_time,
        presolve_bt,
        presolve_timeout,
        presolve_maxiter,
        presolve_bt_width_tol,
        presolve_bt_output_tol,
        presolve_bt_algo,
        presolve_bt_relax,
        presolve_bt_mip_timeout,
        presolve_bp,
        user_parameters,
        int_enable,
        int_cumulative_disc,
        int_fully_disc)
    end

# Create Alpine's nonlinear model -- can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)

    
    # Translate options into old nonlinearmodel.jl fields
    colorful_pod = s.colorful_pod

    loglevel = s.loglevel
    timeout = s.timeout
    maxiter = s.maxiter
    relgap = s.relgap
    gapref = s.gapref
    absgap = s.absgap
    tol = s.tol
    largebound = s.largebound

    recognize_convex = s.recognize_convex
    bilinear_mccormick = s.bilinear_mccormick
    bilinear_convexhull = s.bilinear_convexhull
    monomial_convexhull = s.monomial_convexhull

    method_convexification = s.method_convexification
    method_partition_injection = s.method_partition_injection
    term_patterns = s.term_patterns
    constr_patterns = s.constr_patterns

    nlp_solver = s.nlp_solver
    minlp_solver = s.minlp_solver
    mip_solver = s.mip_solver

    disc_var_pick = s.disc_var_pick
    disc_ratio = s.disc_ratio
    disc_uniform_rate = s.disc_uniform_rate
    disc_add_partition_method = s.disc_add_partition_method
    disc_divert_chunks = s.disc_divert_chunks
    disc_abs_width_tol = s.disc_abs_width_tol
    disc_rel_width_tol = s.disc_rel_width_tol
    disc_consecutive_forbid = s.disc_consecutive_forbid
    disc_ratio_branch = s.disc_ratio_branch

    convhull_formulation = s.convhull_formulation
    convhull_ebd = s.convhull_ebd
    convhull_ebd_encode = s.convhull_ebd_encode
    convhull_ebd_ibs = s.convhull_ebd_ibs
    convhull_ebd_link = s.convhull_ebd_link
    convhull_warmstart = s.convhull_warmstart
    convhull_no_good_cuts = s.convhull_no_good_cuts

    presolve_track_time = s.presolve_track_time
    presolve_bt = s.presolve_bt
    presolve_timeout = s.presolve_timeout
    presolve_maxiter = s.presolve_maxiter
    presolve_bt_width_tol = s.presolve_bt_width_tol
    presolve_bt_output_tol = s.presolve_bt_output_tol
    presolve_bt_algo = s.presolve_bt_algo
    presolve_bt_relax = s.presolve_bt_relax
    presolve_bt_mip_timeout = s.presolve_bt_mip_timeout

    presolve_bp = s.presolve_bp

    user_parameters = s.user_parameters
    int_enable = s.int_enable
    int_cumulative_disc = s.int_cumulative_disc
    int_fully_disc = s.int_fully_disc

    return PODNonlinearModel(colorful_pod,
                            loglevel, timeout, maxiter, relgap, gapref, absgap, tol, largebound,
                            nlp_solver,
                            minlp_solver,
                            mip_solver,
                            recognize_convex,
                            bilinear_mccormick,
                            bilinear_convexhull,
                            monomial_convexhull,
                            method_convexification,
                            method_partition_injection,
                            term_patterns,
                            constr_patterns,
                            disc_var_pick,
                            disc_ratio,
                            disc_uniform_rate,
                            disc_add_partition_method,
                            disc_divert_chunks,
                            disc_abs_width_tol,
                            disc_rel_width_tol,
                            disc_consecutive_forbid,
                            disc_ratio_branch,
                            convhull_formulation,
                            convhull_ebd,
                            convhull_ebd_encode,
                            convhull_ebd_ibs,
                            convhull_ebd_link,
                            convhull_warmstart,
                            convhull_no_good_cuts,
                            presolve_track_time,
                            presolve_bt,
                            presolve_timeout,
                            presolve_maxiter,
                            presolve_bt_width_tol,
                            presolve_bt_output_tol,
                            presolve_bt_algo,
                            presolve_bt_relax,
                            presolve_bt_mip_timeout,
                            presolve_bp,
                            user_parameters,
                            int_enable,
                            int_cumulative_disc,
                            int_fully_disc)
end

function MathProgBase.loadproblem!(m::PODNonlinearModel,
                                   num_var::Int,
                                   num_constr::Int,
                                   l_var::Vector{Float64},
                                   u_var::Vector{Float64},
                                   l_constr::Vector{Float64},
                                   u_constr::Vector{Float64},
                                   sense::Symbol,
                                   d::MathProgBase.AbstractNLPEvaluator)

    # Populating PODNonlinearModel (invoked by JuMP.build(m))
    m.num_var_orig = num_var
    m.num_constr_orig = num_constr
    m.l_var_orig = l_var
    m.u_var_orig = u_var
    m.l_constr_orig = l_constr
    m.u_constr_orig = u_constr
    m.sense_orig = sense
    if m.sense_orig == :Max
        m.best_obj = -Inf
        m.best_bound = Inf
    else
        m.best_obj = Inf
        m.best_bound = -Inf
    end
    m.d_orig = d

    # Initialize NLP interface
    interface_init_nonlinear_data(m.d_orig)

    # Collect objective & constraint expressions
    m.obj_expr_orig = interface_get_obj_expr(m.d_orig)
    for i in 1:m.num_constr_orig
        push!(m.constr_expr_orig, interface_get_constr_expr(m.d_orig, i))
    end

    # Collect original variable type and build dynamic variable type space
    m.var_type_orig = [getcategory(Variable(d.m, i)) for i in 1:m.num_var_orig]
    m.var_type = copy(m.var_type_orig)
    m.int_vars = [i for i in 1:m.num_var_orig if m.var_type[i] == :Int]
    m.bin_vars = [i for i in 1:m.num_var_orig if m.var_type[i] == :Bin]

    if !isempty(m.int_vars) || !isempty(m.bin_vars) 
        (m.minlp_solver == empty_solver) && (error("Problem is a MINLP and no MINLP local solver given; use minlp_solver to specify a MINLP local solver"))
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

    # populate data to create the bounding model
    recategorize_var(m)             # Initial round of variable re-categorization

    :Int in m.var_type_orig && @warn "POD's support for integer variables is experimental"
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

    ## Citation
    println("***********************************************************************")
    println(" This package contains Alpine.jl, a global solver for nonconvex MINLPs")
    println("       If you find it useful, please cite the following paper: ")
    println("     Journal of Global Optimization, 2019, https://goo.gl/89zrDf")
    println("***********************************************************************")

    return
end
