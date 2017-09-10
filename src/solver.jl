export PODSolver

# Dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

type PODSolver <: MathProgBase.AbstractMathProgSolver
    dev_debug::Bool
    dev_test::Bool
    colorful_pod::Any

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
    monomial_convexhull::Bool

    method_convexification::Array{Function}
    expr_patterns::Array{Function}

    discretization_var_pick_algo::Any
    discretization_ratio::Any
    discretization_uniform_rate::Int
    discretization_add_partition_method::Any
    discretization_abs_width_tol::Float64
    discretization_rel_width_tol::Float64
    discretization_consecutive_forbid::Int
    disc_ratio_branch::Bool
    disc_ratio_branch_timeout::Float64
    disc_ratio_branch_focus::Any

    convexhull_sweep_limit::Int
    convhull_formulation_sos2::Bool
    convhull_formulation_sos2aux::Bool
    convhull_formulation_facet::Bool
    convhull_formulation_minib::Bool

    presolve_track_time::Bool
    presolve_bound_tightening::Bool
    presolve_maxiter::Int
    presolve_bt_width_tol::Float64
    presolve_bt_output_tol::Float64
    presolve_bound_tightening_algo::Any
    presolve_mip_relaxation::Bool
    presolve_mip_timelimit::Float64

    bound_basic_propagation::Bool

    user_parameters::Dict

    # other options to be added later on
end

function PODSolver(;
    dev_debug = false,
    dev_test = false,
    colorful_pod = false,

    log_level = 1,
    timeout = Inf,
    maxiter = 99,
    rel_gap = 1e-4,
    tol = 1e-6,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),

    bilinear_mccormick = false,      # by default, deal with bilinear terms using mccormick since the convhull is under tests
    bilinear_convexhull = true,
    monomial_convexhull = true,

    method_convexification = Array{Function}(0),
    expr_patterns = Array{Function}(0),

    discretization_var_pick_algo = 0,           # By default pick all variables
    discretization_ratio = 4,
    discretization_uniform_rate = 2,
    discretization_add_partition_method = "adaptive",
    discretization_abs_width_tol = 1,
    discretization_rel_width_tol = 1e-6,
    discretization_consecutive_forbid = 0,

    disc_ratio_branch=false,
    disc_ratio_branch_timeout=1,
    disc_ratio_branch_focus="gap",

    convexhull_sweep_limit = 1,
    convhull_formulation_sos2 = true,
    convhull_formulation_sos2aux = false,
    convhull_formulation_facet = false,
    convhull_formulation_minib = false,

    presolve_track_time = false,
    presolve_bound_tightening = false,
    presolve_maxiter = 9999,
    presolve_bt_width_tol = 1e-3,
    presolve_bt_output_tol = 1e-5,
    presolve_bound_tightening_algo = 1,
    presolve_mip_relaxation = false,
    presolve_mip_timelimit = Inf,

    bound_basic_propagation = false,
    user_parameters = Dict(),
    kwargs...
    )

    unsupport_opts = Dict(kwargs)
    !isempty(keys(unsupport_opts)) && warn("Detected unsupported/experimental arguments = $(keys(unsupported_opts))")

    if nlp_local_solver == UnsetSolver()
        error("No NLP local solver specified (set nlp_local_solver)\n")
    end

    if mip_solver == UnsetSolver()
        error("NO MIP solver specififed (set mip_solver)\n")
    end

    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(dev_debug, dev_test, colorful_pod,
        log_level, timeout, maxiter, rel_gap, tol,
        deepcopy(nlp_local_solver),
        deepcopy(minlp_local_solver),
        deepcopy(mip_solver),
        bilinear_mccormick,
        bilinear_convexhull,
        monomial_convexhull,
        method_convexification,
        expr_patterns,
        discretization_var_pick_algo,
        discretization_ratio,
        discretization_uniform_rate,
        discretization_add_partition_method,
        discretization_abs_width_tol,
        discretization_rel_width_tol,
        discretization_consecutive_forbid,
        disc_ratio_branch,
        disc_ratio_branch_timeout,
        disc_ratio_branch_focus,
        convexhull_sweep_limit,
        convhull_formulation_sos2,
        convhull_formulation_sos2aux,
        convhull_formulation_facet,
        convhull_formulation_minib,
        presolve_track_time,
        presolve_bound_tightening,
        presolve_maxiter,
        presolve_bt_width_tol,
        presolve_bt_output_tol,
        presolve_bound_tightening_algo,
        presolve_mip_relaxation,
        presolve_mip_timelimit,
        bound_basic_propagation,
        user_parameters)
    end

# Create POD nonlinear model: can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)
    if !applicable(MathProgBase.NonlinearModel, s.nlp_local_solver)
        error("NLP local solver $(s.nlp_local_solver) specified is not a NLP solver recognized by POD\n")
    end

    # Translate options into old nonlinearmodel.jl fields
    dev_test = s.dev_test
    dev_debug = s.dev_debug
    colorful_pod = s.colorful_pod

    log_level = s.log_level
    timeout = s.timeout
    maxiter = s.maxiter
    rel_gap = s.rel_gap
    tol = s.tol

    bilinear_mccormick = s.bilinear_mccormick
    bilinear_convexhull = s.bilinear_convexhull
    monomial_convexhull = s.monomial_convexhull

    method_convexification = s.method_convexification
    expr_patterns = s.expr_patterns

    nlp_local_solver = s.nlp_local_solver
    minlp_local_solver = s.minlp_local_solver
    mip_solver = s.mip_solver

    discretization_var_pick_algo = s.discretization_var_pick_algo
    discretization_ratio = s.discretization_ratio
    discretization_uniform_rate = s.discretization_uniform_rate
    discretization_add_partition_method = s.discretization_add_partition_method
    discretization_abs_width_tol = s.discretization_abs_width_tol
    discretization_rel_width_tol = s.discretization_rel_width_tol
    discretization_consecutive_forbid = s.discretization_consecutive_forbid

    disc_ratio_branch = s.disc_ratio_branch
    disc_ratio_branch_timeout = s.disc_ratio_branch_timeout
    disc_ratio_branch_focus = s.disc_ratio_branch_focus

    convexhull_sweep_limit = s.convexhull_sweep_limit
    convhull_formulation_sos2 = s.convhull_formulation_sos2
    convhull_formulation_sos2aux = s.convhull_formulation_sos2aux
    convhull_formulation_facet = s.convhull_formulation_facet
    convhull_formulation_minib = s.convhull_formulation_minib

    presolve_track_time = s.presolve_track_time
    presolve_bound_tightening = s.presolve_bound_tightening
    presolve_maxiter = s.presolve_maxiter
    presolve_bt_width_tol = s.presolve_bt_width_tol
    presolve_bt_output_tol = s.presolve_bt_output_tol
    presolve_bound_tightening_algo = s.presolve_bound_tightening_algo
    presolve_mip_relaxation = s.presolve_mip_relaxation
    presolve_mip_timelimit = s.presolve_mip_timelimit

    bound_basic_propagation = s.bound_basic_propagation

    user_parameters = s.user_parameters

    return PODNonlinearModel(dev_debug, dev_test, colorful_pod,
                            log_level, timeout, maxiter, rel_gap, tol,
                            nlp_local_solver,
                            minlp_local_solver,
                            mip_solver,
                            bilinear_mccormick,
                            bilinear_convexhull,
                            monomial_convexhull,
                            method_convexification,
                            expr_patterns,
                            discretization_var_pick_algo,
                            discretization_ratio,
                            discretization_uniform_rate,
                            discretization_add_partition_method,
                            discretization_abs_width_tol,
                            discretization_rel_width_tol,
                            discretization_consecutive_forbid,
                            disc_ratio_branch,
                            disc_ratio_branch_timeout,
                            disc_ratio_branch_focus,
                            convexhull_sweep_limit,
                            convhull_formulation_sos2,
                            convhull_formulation_sos2aux,
                            convhull_formulation_facet,
                            convhull_formulation_minib,
                            presolve_track_time,
                            presolve_bound_tightening,
                            presolve_maxiter,
                            presolve_bt_width_tol,
                            presolve_bt_output_tol,
                            presolve_bound_tightening_algo,
                            presolve_mip_relaxation,
                            presolve_mip_timelimit,
                            bound_basic_propagation,
                            user_parameters)
end
