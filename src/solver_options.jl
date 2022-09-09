# Optimizer struct and user-defined options (with default settings) to tune Alpine solver

mutable struct OptimizerOptions
    # Parameters for tuning Alpine

    # Basic solver parameters
    log_level::Int                                              # Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info
    time_limit::Float64                                         # Time limit for algorithm (in seconds)
    max_iter::Int                                               # Maximum bound on number of iterations for the partitioning algorithm
    rel_gap::Float64                                            # Relative optimality gap for termination condition
    abs_gap::Float64                                            # Absolute optimality gap for termination condition
    tol::Float64                                                # Numerical tol used in algorithms
    large_bound::Float64                                        # Large bounds for problems with unbounded variables

    # All the solver options
    nlp_solver::Union{MOI.OptimizerWithAttributes,Nothing}      # Local continuous NLP solver for solving NLPs at each iteration
    minlp_solver::Union{MOI.OptimizerWithAttributes,Nothing}    # Local MINLP solver for solving MINLPs at each iteration
    mip_solver::Union{MOI.OptimizerWithAttributes,Nothing}      # MIP solver for successive lower bound solves

    # Convexification methods
    recognize_convex::Bool                                      # Recognize convex expressions in parsing objective functions and constraints

    # Expression-based user-inputs
    # method_convexification :: Array{Function}                 # Array of functions that user can choose to convexify specific non-linear terms : no over-ride privilege

    # Parameters used in Alpine's MIP-based partitioning algorithm
    apply_partitioning::Bool                                    # Apply the partitioning algorithm only if thhis true, else terminate after presolve
    disc_var_pick::Any                                          # Algorithm for choosing the variables to discretize: 1 for minimum vertex cover, 0 for all variables
    partition_scaling_factor::Int                               # Partition scaling parameter, which is critical for convergence (using a fixed value for now, later switch to a function)
    disc_uniform_rate::Int                                      # Discretization rate parameter when using uniform partitions
    disc_add_partition_method::Any                              # Additional methods to add discretization
    disc_divert_chunks::Int                                     # How many uniform partitions to construct
    disc_abs_width_tol::Float64                                 # Absolute tolerance used when setting up a partition
    disc_rel_width_tol::Float64                                 # Relative width tolerance when setting up a partition
    disc_consecutive_forbid::Int                                # Prevent bounding model to add partitions consecutively in the same region when bounds do not improve
    partition_scaling_factor_branch::Bool                       # Branching tests for picking fixed discretization ratio

    # MIP formulation parameters
    convhull_formulation::String                                # MIP Formulation for the relaxation
    convhull_ebd::Bool                                          # Enable embedding formulation
    convhull_ebd_encode::Any                                    # Encoding method used for convhull_ebd
    convhull_ebd_ibs::Bool                                      # Enable independent branching scheme
    convhull_ebd_link::Bool                                     # Linking constraints between x and α, type 1 uses hierarchical and type 2 uses big-m
    convhull_warmstart::Bool                                    # Warm start the bounding MIP
    linking_constraints::Bool                                   # Multilinear linking constraint feature
    linking_constraints_degree_limit::Int64                     # Degree of shared multilinear terms up to which the linking constraints are built and added

    # Presolve and bound-tightening parameters
    presolve_track_time::Bool                                   # Account presolve time for total time usage
    presolve_bt::Bool                                           # Perform bound tightening procedure before the main algorithm (default: true)
    presolve_bt_time_limit::Float64                             # Time limit for presolving (seconds)
    presolve_bt_max_iter::Int                                   # Maximum iterations allowed to perform presolve
    presolve_bt_width_tol::Float64                              # Width tolerance for bound-tightening
    presolve_bt_improv_tol::Float64                             # Improvement tolerance for average reduction of variable ranges
    presolve_bt_bound_tol::Float64                              # Variable bounds truncation precicision tol
    presolve_bt_obj_bound_tol::Float64                          # Objective upper bound truncation tol
    presolve_bt_algo::Any                                       # Method used for bound tightening procedures, can either be an index of default methods or functional inputs
    presolve_bt_relax_integrality::Bool                         # Relax the MIP solved in built-in relaxation scheme for time performance
    presolve_bt_mip_time_limit::Float64                         # Time limit for a single MIP solved in the built-in bound tightening algorithm (with partitions)
    use_start_as_local_solution::Bool                           # Use starting value as local solution for presolve without invoking a local solver

    # Domain Reduction
    presolve_bp::Bool                                           # Conduct basic bound propagation
end

function get_default_options()
    log_level = 1
    time_limit = 1E6
    max_iter = 99
    rel_gap = 1e-4
    abs_gap = 1e-6
    tol = 1e-6
    large_bound = 1E6

    nlp_solver = nothing
    minlp_solver = nothing
    mip_solver = nothing

    recognize_convex = true

    apply_partitioning = true
    disc_var_pick = 2                # By default, uses the 15-variable selective rule
    partition_scaling_factor = 10
    disc_uniform_rate = 2
    disc_add_partition_method = "adaptive"
    disc_divert_chunks = 5
    disc_abs_width_tol = 1e-4
    disc_rel_width_tol = 1e-6
    disc_consecutive_forbid = false
    partition_scaling_factor_branch = false

    convhull_formulation = "sos2"
    convhull_ebd = false
    convhull_ebd_encode = "default"
    convhull_ebd_ibs = false
    convhull_ebd_link = false
    convhull_warmstart = true
    linking_constraints = true
    linking_constraints_degree_limit = 4 # check up to quadrilinear sharing across terms

    presolve_track_time = true
    presolve_bt = true
    presolve_bt_time_limit = 900
    presolve_bt_max_iter = 25
    presolve_bt_width_tol = 1e-2
    presolve_bt_improv_tol = 1e-3
    presolve_bt_bound_tol = 1e-4
    presolve_bt_obj_bound_tol = 1e-2
    presolve_bt_algo = 1
    presolve_bt_relax_integrality = false
    presolve_bt_mip_time_limit = Inf
    use_start_as_local_solution = false
    presolve_bp = false

    return OptimizerOptions(
        log_level,
        time_limit,
        max_iter,
        rel_gap,
        abs_gap,
        tol,
        large_bound,
        nlp_solver,
        minlp_solver,
        mip_solver,
        recognize_convex,
        apply_partitioning,
        disc_var_pick,
        partition_scaling_factor,
        disc_uniform_rate,
        disc_add_partition_method,
        disc_divert_chunks,
        disc_abs_width_tol,
        disc_rel_width_tol,
        disc_consecutive_forbid,
        partition_scaling_factor_branch,
        convhull_formulation,
        convhull_ebd,
        convhull_ebd_encode,
        convhull_ebd_ibs,
        convhull_ebd_link,
        convhull_warmstart,
        linking_constraints,
        linking_constraints_degree_limit,
        presolve_track_time,
        presolve_bt,
        presolve_bt_time_limit,
        presolve_bt_max_iter,
        presolve_bt_width_tol,
        presolve_bt_improv_tol,
        presolve_bt_bound_tol,
        presolve_bt_obj_bound_tol,
        presolve_bt_algo,
        presolve_bt_relax_integrality,
        presolve_bt_mip_time_limit,
        use_start_as_local_solution,
        presolve_bp,
    )
end
