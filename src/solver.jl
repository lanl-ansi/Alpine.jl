# Parameters for tuning Alpine, a global solver for non-convex MINLPs 

mutable struct OptimizerOptions

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
    presolve_bt_relax_integrality::Bool                         # Relax the MIP solved in built-in relaxation scheme for time performance
    presolve_bt_mip_time_limit::Float64                         # Time limit for a single MIP solved in the built-in bound tightening algorithm (with partitions)

    # Domain Reduction
    presolve_bp::Bool                                           # Conduct basic bound propagation
    user_parameters::Dict                                       # [INACTIVE] Additional parameters used for user-defined functional inputs

    # Features for Integer Problems (NOTE: no support for int-lin problems)
    int_enable::Bool                                            # Convert integer problem into binary problem
    int_cumulative_disc::Bool                                   # Cumulatively involve integer variables for discretization
    
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

    disc_var_pick = 2                      # Check pick_disc_vars function for more options
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
    presolve_bt_max_iter = 5
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