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
    max_iter::Int
    rel_gap::Float64
    tol::Float64

    nlp_local_solver::MathProgBase.AbstractMathProgSolver
    minlp_local_solver::MathProgBase.AbstractMathProgSolver
    mip_solver::MathProgBase.AbstractMathProgSolver

    recognize_convex::Bool
    bilinear_mccormick::Bool
    bilinear_convexhull::Bool
    monomial_convexhull::Bool

    method_convexification::Array{Function}
    term_patterns::Array{Function}
    constr_patterns::Array{Function}

    disc_var_pick_algo::Any
    disc_ratio::Any
    disc_uniform_rate::Int
    disc_add_partition_method::Any
    disc_abs_width_tol::Float64
    disc_rel_width_tol::Float64
    disc_consecutive_forbid::Int
    disc_ratio_branch::Bool

    convexhull_sweep_limit::Int
    convhull_formulation_sos2::Bool
    convhull_formulation_sos2aux::Bool
    convhull_formulation_facet::Bool
    convhull_formulation_minib::Bool

    presolve_track_time::Bool
    presolve_bound_tightening::Bool
    presolve_max_iter::Int
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
    max_iter = 99,
    rel_gap = 1e-4,
    tol = 1e-6,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),

    recognize_convex = true,
    bilinear_mccormick = false,
    bilinear_convexhull = true,
    monomial_convexhull = true,

    method_convexification = Array{Function}(0),
    term_patterns = Array{Function}(0),
    constr_patterns = Array{Function}(0),

    disc_var_pick_algo = 0,           # By default pick all variables
    disc_ratio = 4,
    disc_uniform_rate = 2,
    disc_add_partition_method = "adaptive",
    disc_abs_width_tol = 1e-4,
    disc_rel_width_tol = 1e-6,
    disc_consecutive_forbid = 0,
    disc_ratio_branch=false,

    convexhull_sweep_limit = 1,
    convhull_formulation_sos2 = true,
    convhull_formulation_sos2aux = false,
    convhull_formulation_facet = false,
    convhull_formulation_minib = false,

    presolve_track_time = true,
    presolve_bound_tightening = false,
    presolve_max_iter = 9999,
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
    !isempty(keys(unsupport_opts)) && warn("Detected unsupported/experimental arguments = $(keys(unsupport_opts))")

    if nlp_local_solver == UnsetSolver()
        error("No NLP local solver specified (set nlp_local_solver)\n")
    end

    if mip_solver == UnsetSolver()
        error("NO MIP solver specififed (set mip_solver)\n")
    end

    # Code Conversion
    (discretization_var_pick_algo in ["select_all_nlvar", "all", "max"]) && (discretization_var_pick_algo = 0)
    (discretization_var_pick_algo in ["min_vertex_cover","min"]) && (discretization_var_pick_algo = 1)
    (discretization_var_pick_algo == "selective") && (discretization_var_pick_algo = 2)
    (discretization_var_pick_algo == "dynamic") && (discretization_var_pick_algo = 3)

    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(dev_debug, dev_test, colorful_pod,
        log_level, timeout, max_iter, rel_gap, tol,
        deepcopy(nlp_local_solver),
        deepcopy(minlp_local_solver),
        deepcopy(mip_solver),
        recognize_convex,
        bilinear_mccormick,
        bilinear_convexhull,
        monomial_convexhull,
        method_convexification,
        term_patterns,
        constr_patterns,
        disc_var_pick_algo,
        disc_ratio,
        disc_uniform_rate,
        disc_add_partition_method,
        disc_abs_width_tol,
        disc_rel_width_tol,
        disc_consecutive_forbid,
        disc_ratio_branch,
        convexhull_sweep_limit,
        convhull_formulation_sos2,
        convhull_formulation_sos2aux,
        convhull_formulation_facet,
        convhull_formulation_minib,
        presolve_track_time,
        presolve_bound_tightening,
        presolve_max_iter,
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
    max_iter = s.max_iter
    rel_gap = s.rel_gap
    tol = s.tol

    recognize_convex = s.recognize_convex
    bilinear_mccormick = s.bilinear_mccormick
    bilinear_convexhull = s.bilinear_convexhull
    monomial_convexhull = s.monomial_convexhull

    method_convexification = s.method_convexification
    term_patterns = s.term_patterns
    constr_patterns = s.constr_patterns

    nlp_local_solver = s.nlp_local_solver
    minlp_local_solver = s.minlp_local_solver
    mip_solver = s.mip_solver

    disc_var_pick_algo = s.disc_var_pick_algo
    disc_ratio = s.disc_ratio
    disc_uniform_rate = s.disc_uniform_rate
    disc_add_partition_method = s.disc_add_partition_method
    disc_abs_width_tol = s.disc_abs_width_tol
    disc_rel_width_tol = s.disc_rel_width_tol
    disc_consecutive_forbid = s.disc_consecutive_forbid
    disc_ratio_branch = s.disc_ratio_branch

    convexhull_sweep_limit = s.convexhull_sweep_limit
    convhull_formulation_sos2 = s.convhull_formulation_sos2
    convhull_formulation_sos2aux = s.convhull_formulation_sos2aux
    convhull_formulation_facet = s.convhull_formulation_facet
    convhull_formulation_minib = s.convhull_formulation_minib

    presolve_track_time = s.presolve_track_time
    presolve_bound_tightening = s.presolve_bound_tightening
    presolve_max_iter = s.presolve_max_iter
    presolve_bt_width_tol = s.presolve_bt_width_tol
    presolve_bt_output_tol = s.presolve_bt_output_tol
    presolve_bound_tightening_algo = s.presolve_bound_tightening_algo
    presolve_mip_relaxation = s.presolve_mip_relaxation
    presolve_mip_timelimit = s.presolve_mip_timelimit

    bound_basic_propagation = s.bound_basic_propagation

    user_parameters = s.user_parameters

    return PODNonlinearModel(dev_debug, dev_test, colorful_pod,
                            log_level, timeout, max_iter, rel_gap, tol,
                            nlp_local_solver,
                            minlp_local_solver,
                            mip_solver,
                            recognize_convex,
                            bilinear_mccormick,
                            bilinear_convexhull,
                            monomial_convexhull,
                            method_convexification,
                            term_patterns,
                            constr_patterns,
                            disc_var_pick_algo,
                            disc_ratio,
                            disc_uniform_rate,
                            disc_add_partition_method,
                            disc_abs_width_tol,
                            disc_rel_width_tol,
                            disc_consecutive_forbid,
                            disc_ratio_branch,
                            convexhull_sweep_limit,
                            convhull_formulation_sos2,
                            convhull_formulation_sos2aux,
                            convhull_formulation_facet,
                            convhull_formulation_minib,
                            presolve_track_time,
                            presolve_bound_tightening,
                            presolve_max_iter,
                            presolve_bt_width_tol,
                            presolve_bt_output_tol,
                            presolve_bound_tightening_algo,
                            presolve_mip_relaxation,
                            presolve_mip_timelimit,
                            bound_basic_propagation,
                            user_parameters)
end

# MathProgBase.LinearQuadraticModel(s::PODSolver) = MathProgBase.NonlinearModel(s::PODSolver)

function fetch_mip_solver_identifier(m::PODNonlinearModel)

    if string(m.mip_solver)[1:6] == "Gurobi"
        m.mip_solver_identifier = "Gurobi"
    elseif string(m.mip_solver)[1:5] == "CPLEX"
        m.mip_solver_identifier = "CPLEX"
    elseif string(m.mip_solver)[1:3] == "Cbc"
        m.mip_solver_identifier = "Cbc"
    elseif string(m.mip_solver)[1:4] == "GLPK"
        m.mip_solver_identifier = "GLPK"
    elseif string(m.mip_solver)[1:8] == "Pajarito"
        m.mip_solver_identifier = "Pajarito"
    else
        error("Unsupported mip solver name. Using blank")
    end

    return
end

function fetch_nlp_solver_identifier(m::PODNonlinearModel)

    if string(m.nlp_local_solver)[1:5] == "Ipopt"
        m.nlp_local_solver_identifier = "Ipopt"
    elseif string(m.nlp_local_solver)[1:6] == "AmplNL"
        m.nlp_local_solver_identifier = "Bonmin"
    elseif string(m.nlp_local_solver)[1:6] == "Knitro"
        m.nlp_local_solver_identifier = "Knitro"
    elseif string(m.nlp_local_solver)[1:8] == "Pajarito"
        m.nlp_local_solver_identifier = "Pajarito"
    elseif string(m.nlp_local_solver)[1:5] == "NLopt"
        m.nlp_local_solver_identifier = "NLopt"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

function fetch_minlp_solver_identifier(m::PODNonlinearModel)

    (m.minlp_local_solver == UnsetSolver()) && return
    if string(m.minlp_local_solver)[1:6] == "AmplNL"
        m.minlp_local_solver_identifier = "Bonmin"
    elseif string(m.minlp_local_solver)[1:6] == "Knitro"
        m.minlp_local_solver_identifier = "Knitro"
    elseif string(m.minlp_local_solver)[1:8] == "Pajarito"
        m.minlp_local_solver_identifier = "Pajarito"
    elseif string(m.minlp_local_solver)[1:5] == "NLopt"
        m.minlp_local_solver_identifier = "NLopt"
    elseif string(m.minlp_local_solver)[1:35] == "CoinOptServices.OsilSolver(\"bonmin\""
        m.minlp_local_solver_identifier = "Bonmin"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_mip_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.mip_solver_identifier == "CPLEX"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.mip_solver_identifier == "Gurobi"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:TimeLimit,m.timeout)
    elseif m.mip_solver_identifier == "Cbc"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:seconds,m.timeout)
    elseif m.mip_solver_identifier == "GLPK"
        insert_timeleft_symbol(m.mip_solver.opt)
    elseif m.mip_solver_identifier == "Pajarito"
        (timelimit < Inf) && (m.mip_solver.timeout = timelimit)
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_nlp_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.nlp_local_solver_identifier == "Ipopt"
        insert_timeleft_symbol(m.nlp_local_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.nlp_local_solver_identifier == "Pajarito"
        (timelimit < Inf) && (m.nlp_local_solver.timeout = timelimit)
    elseif m.nlp_local_solver_identifier == "AmplNL"
        insert_timeleft_symbol(m.nlp_local_solver.options,timelimit,:seconds,m.timeout, options_string_type=2)
    elseif m.nlp_local_solver_identifier == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    elseif m.nlp_local_solver_identifier == "NLopt"
        m.nlp_local_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_minlp_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.minlp_local_solver_identifier == "Pajarito"
        (timelimit < Inf) && (m.minlp_local_solver.timeout = timelimit)
    elseif m.minlp_local_solver_identifier == "AmplNL"
        insert_timeleft_symbol(m.minlp_local_solver.options,timelimit,:seconds,m.timeout,options_string_type=2)
    elseif m.minlp_local_solver_identifier == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    elseif m.minlp_local_solver_identifier == "NLopt"
        m.minlp_local_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end
