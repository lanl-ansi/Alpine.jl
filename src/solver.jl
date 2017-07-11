export PODSolver

# Dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

type PODSolver <: MathProgBase.AbstractMathProgSolver
    dev_debug::Bool
    dev_test::Bool

    log_level::Int
    timeout::Float64
    maxiter::Int
    rel_gap::Float64
    tol::Float64

    nlp_local_solver::MathProgBase.AbstractMathProgSolver
    minlp_local_solver::MathProgBase.AbstractMathProgSolver
    mip_solver::MathProgBase.AbstractMathProgSolver

    bilinear_mccormick::Bool
    bilinear_convexhull::Bool

    method_convexification::Array{Function}
    expr_patterns::Array{Function}

    discretization_var_pick_algo::Any
    discretization_ratio::Any
    discretization_add_partition_method::Any

    presolve_track_time::Bool
    presolve_bound_tightening::Bool
    presolve_maxiter::Int
    presolve_bt_width_tol::Float64
    presolve_bt_output_tol::Float64
    presolve_bound_tightening_algo::Any
    presolve_mip_relaxation::Bool
    presolve_mip_timelimit::Float64

    # other options to be added later on
end

function PODSolver(;
    dev_debug = false,
    dev_test = false,

    log_level = 1,
    timeout = Inf,
    maxiter = 99,
    rel_gap = 1e-4,
    tol = 1e-6,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),

    bilinear_mccormick = true,      # by default, deal with bilinear terms using mccormick
    bilinear_convexhull = false,

    method_convexification = Array{Function}(0),
    expr_patterns = Array{Function}(0),

    discretization_var_pick_algo = 0,           # By default pick all variables
    discretization_ratio = 4,
    discretization_add_partition_method = nothing, # Not ready for implementation

    presolve_track_time = false,
    presolve_bound_tightening = false,
    presolve_maxiter = 9999,
    presolve_bt_width_tol = 1e-3,
    presolve_bt_output_tol = 1e-5,
    presolve_bound_tightening_algo = 1,
    presolve_mip_relaxation = false,
    presolve_mip_timelimit = Inf,
    )

    if nlp_local_solver == UnsetSolver()
        error("No NLP local solver specified (set nlp_local_solver)\n")
    end

    if mip_solver == UnsetSolver()
        error("NO MIP solver specififed (set mip_solver)\n")
    end

    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(dev_debug, dev_test,
        log_level, timeout, maxiter, rel_gap, tol,
        deepcopy(nlp_local_solver),
        deepcopy(minlp_local_solver),
        deepcopy(mip_solver),
        bilinear_mccormick,
        bilinear_convexhull,
        method_convexification,
        expr_patterns,
        discretization_var_pick_algo,
        discretization_ratio,
        discretization_add_partition_method,
        presolve_track_time,
        presolve_bound_tightening,
        presolve_maxiter,
        presolve_bt_width_tol,
        presolve_bt_output_tol,
        presolve_bound_tightening_algo,
        presolve_mip_relaxation,presolve_mip_timelimit)
    end

# Create POD nonlinear model: can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)
    if !applicable(MathProgBase.NonlinearModel, s.nlp_local_solver)
        error("NLP local solver $(s.nlp_local_solver) specified is not a NLP solver recognized by POD\n")
    end

    # Translate options into old nonlinearmodel.jl fields
    dev_test = s.dev_test
    dev_debug = s.dev_debug

    log_level = s.log_level
    timeout = s.timeout
    maxiter = s.maxiter
    rel_gap = s.rel_gap
    tol = s.tol

    bilinear_mccormick = s.bilinear_mccormick
    bilinear_convexhull = s.bilinear_convexhull

    method_convexification = s.method_convexification
    expr_patterns = s.expr_patterns

    nlp_local_solver = s.nlp_local_solver
    minlp_local_solver = s.minlp_local_solver
    mip_solver = s.mip_solver

    discretization_var_pick_algo = s.discretization_var_pick_algo
    discretization_ratio = s.discretization_ratio
    discretization_add_partition_method = s.discretization_add_partition_method

    presolve_track_time = s.presolve_track_time
    presolve_bound_tightening = s.presolve_bound_tightening
    presolve_maxiter = s.presolve_maxiter
    presolve_bt_width_tol = s.presolve_bt_width_tol
    presolve_bt_output_tol = s.presolve_bt_output_tol
    presolve_bound_tightening_algo = s.presolve_bound_tightening_algo
    presolve_mip_relaxation = s.presolve_mip_relaxation
    presolve_mip_timelimit = s.presolve_mip_timelimit

    return PODNonlinearModel(dev_debug, dev_test,
                            log_level, timeout, maxiter, rel_gap, tol,
                            nlp_local_solver,
                            minlp_local_solver,
                            mip_solver,
                            bilinear_mccormick,
                            bilinear_convexhull,
                            method_convexification,
                            expr_patterns,
                            discretization_var_pick_algo,
                            discretization_ratio,
                            discretization_add_partition_method,
                            presolve_track_time,
                            presolve_bound_tightening,
                            presolve_maxiter,
                            presolve_bt_width_tol,
                            presolve_bt_output_tol,
                            presolve_bound_tightening_algo,
                            presolve_mip_relaxation,
                            presolve_mip_timelimit)
end
