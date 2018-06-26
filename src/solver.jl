export PODSolver

# Dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

type PODSolver <: MathProgBase.AbstractMathProgSolver

    colorful_pod::Any

    loglevel::Int
    timeout::Float64
    maxiter::Int
    relgap::Float64
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

    disc_var_pick::Any
    disc_ratio::Any
    disc_uniform_rate::Int
    disc_add_partition_method::Any
    disc_abs_width_tol::Float64
    disc_rel_width_tol::Float64
    disc_consecutive_forbid::Int
    disc_ratio_branch::Bool

    convhull_formulation::String
    convhull_ebd::Bool
    convhull_ebd_encode::Any
    convhull_ebd_ibs::Bool
    convhull_ebd_link::Bool

    presolve_track_time::Bool
    presolve_bt::Bool
    presolve_maxiter::Int
    presolve_bt_width_tol::Float64
    presolve_bt_output_tol::Float64
    presolve_bt_algo::Any
    presolve_bt_relax::Bool
    presolve_bt_mip_timeout::Float64

    presolve_bp::Bool

    user_parameters::Dict

    # other options to be added later on
end

function PODSolver(;

    colorful_pod = false,

    loglevel = 1,
    timeout = Inf,
    maxiter = 99,
    relgap = 1e-4,
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

    disc_var_pick = 2,                     # By default use the 15-variable selective rule
    disc_ratio = 4,
    disc_uniform_rate = 2,
    disc_add_partition_method = "adaptive",
    disc_abs_width_tol = 1e-4,
    disc_rel_width_tol = 1e-6,
    disc_consecutive_forbid = 0,
    disc_ratio_branch=false,

    convhull_formulation = "sos2",
    convhull_ebd = false,
    convhull_ebd_encode = "default",
    convhull_ebd_ibs = false,
    convhull_ebd_link = false,

    presolve_track_time = true,
    presolve_maxiter = 9999,
    presolve_bt = false,
    presolve_bt_width_tol = 1e-3,
    presolve_bt_output_tol = 1e-5,
    presolve_bt_algo = 1,
    presolve_bt_relax = false,
    presolve_bt_mip_timeout = Inf,

    presolve_bp = true,
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
    (disc_var_pick in ["select_all_nlvar", "all", "max"]) && (disc_var_pick = 0)
    (disc_var_pick in ["min_vertex_cover","min"]) && (disc_var_pick = 1)
    (disc_var_pick == "selective") && (disc_var_pick = 2)
    (disc_var_pick == "dynamic") && (disc_var_pick = 3)

    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(colorful_pod,
        loglevel, timeout, maxiter, relgap, tol,
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
        disc_var_pick,
        disc_ratio,
        disc_uniform_rate,
        disc_add_partition_method,
        disc_abs_width_tol,
        disc_rel_width_tol,
        disc_consecutive_forbid,
        disc_ratio_branch,
        convhull_formulation,
        convhull_ebd,
        convhull_ebd_encode,
        convhull_ebd_ibs,
        convhull_ebd_link,
        presolve_track_time,
        presolve_bt,
        presolve_maxiter,
        presolve_bt_width_tol,
        presolve_bt_output_tol,
        presolve_bt_algo,
        presolve_bt_relax,
        presolve_bt_mip_timeout,
        presolve_bp,
        user_parameters)
    end

# Create POD nonlinear model: can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)
    if !applicable(MathProgBase.NonlinearModel, s.nlp_local_solver)
        error("NLP local solver $(s.nlp_local_solver) specified is not a NLP solver recognized by POD\n")
    end

    # Translate options into old nonlinearmodel.jl fields
    colorful_pod = s.colorful_pod

    loglevel = s.loglevel
    timeout = s.timeout
    maxiter = s.maxiter
    relgap = s.relgap
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

    disc_var_pick = s.disc_var_pick
    disc_ratio = s.disc_ratio
    disc_uniform_rate = s.disc_uniform_rate
    disc_add_partition_method = s.disc_add_partition_method
    disc_abs_width_tol = s.disc_abs_width_tol
    disc_rel_width_tol = s.disc_rel_width_tol
    disc_consecutive_forbid = s.disc_consecutive_forbid
    disc_ratio_branch = s.disc_ratio_branch

    convhull_formulation = s.convhull_formulation
    convhull_ebd = s.convhull_ebd
    convhull_ebd_encode = s.convhull_ebd_encode
    convhull_ebd_ibs = s.convhull_ebd_ibs
    convhull_ebd_link = s.convhull_ebd_link

    presolve_track_time = s.presolve_track_time
    presolve_bt = s.presolve_bt
    presolve_maxiter = s.presolve_maxiter
    presolve_bt_width_tol = s.presolve_bt_width_tol
    presolve_bt_output_tol = s.presolve_bt_output_tol
    presolve_bt_algo = s.presolve_bt_algo
    presolve_bt_relax = s.presolve_bt_relax
    presolve_bt_mip_timeout = s.presolve_bt_mip_timeout

    presolve_bp = s.presolve_bp

    user_parameters = s.user_parameters

    return PODNonlinearModel(colorful_pod,
                            loglevel, timeout, maxiter, relgap, tol,
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
                            disc_var_pick,
                            disc_ratio,
                            disc_uniform_rate,
                            disc_add_partition_method,
                            disc_abs_width_tol,
                            disc_rel_width_tol,
                            disc_consecutive_forbid,
                            disc_ratio_branch,
                            convhull_formulation,
                            convhull_ebd,
                            convhull_ebd_encode,
                            convhull_ebd_ibs,
                            convhull_ebd_link,
                            presolve_track_time,
                            presolve_bt,
                            presolve_maxiter,
                            presolve_bt_width_tol,
                            presolve_bt_output_tol,
                            presolve_bt_algo,
                            presolve_bt_relax,
                            presolve_bt_mip_timeout,
                            presolve_bp,
                            user_parameters)
end

# MathProgBase.LinearQuadraticModel(s::PODSolver) = MathProgBase.NonlinearModel(s::PODSolver)

function fetch_mip_solver_identifier(m::PODNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.mip_solver) : solverstring=override

    # Higher-level solvers: that can use sub-solvers
    if contains(solverstring,"Pajarito")
        m.mip_solver_identifier = "Pajarito"
        return
    end

    # Lower level solvers
    if contains(solverstring,"Gurobi")
        m.mip_solver_identifier = "Gurobi"
    elseif contains(solverstring,"CPLEX")
        m.mip_solver_identifier = "CPLEX"
    elseif contains(solverstring,"Cbc")
        m.mip_solver_identifier = "Cbc"
    elseif contains(solverstring,"GLPK")
        m.mip_solver_identifier = "GLPK"
    else
        error("Unsupported mip solver name. Using blank")
    end

    return
end

function fetch_nlp_solver_identifier(m::PODNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.nlp_local_solver) : solverstring=override

    # Higher-level solver
    if contains(solverstring, "Pajarito")
        m.nlp_local_solver_identifier = "Pajarito"
        return
    end

    # Lower-level solver
    if contains(solverstring, "Ipopt")
        m.nlp_local_solver_identifier = "Ipopt"
    elseif contains(solverstring, "AmplNL") && contains(solverstring, "bonmin")
        m.nlp_local_solver_identifier = "Bonmin"
    elseif contains(solverstring, "KNITRO")
        m.nlp_local_solver_identifier = "Knitro"
    elseif contains(solverstring, "NLopt")
        m.nlp_local_solver_identifier = "NLopt"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

function fetch_minlp_solver_identifier(m::PODNonlinearModel;override="")

    (m.minlp_local_solver == UnsetSolver()) && return

    isempty(override) ? solverstring = string(m.minlp_local_solver) : solverstring=override

    # Higher-level solver
    if contains(solverstring, "Pajarito")
        m.minlp_local_solver_identifier = "Pajarito"
        return
    end

    # Lower-level Solver
    if contains(solverstring, "AmplNL") && contains(solverstring, "bonmin")
        m.minlp_local_solver_identifier = "Bonmin"
    elseif contains(solverstring, "KNITRO")
        m.minlp_local_solver_identifier = "Knitro"
    elseif contains(solverstring, "NLopt")
        m.minlp_local_solver_identifier = "NLopt"
    elseif contains(solverstring, "CoinOptServices.OsilSolver(\"bonmin\"")
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
        insert_timeleft_symbol(m.mip_solver.opts, timelimit,:tm_lim,m.timeout)
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
        error("You never tell me anything about knitro. Probably because they have a very short trail length.")
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
