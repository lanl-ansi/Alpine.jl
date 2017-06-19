export PODSolver

# Dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

type PODSolver <: MathProgBase.AbstractMathProgSolver
    log_level::Int
    timeout::Float64
    maxiter::Int
    rel_gap::Float64
    tolerance::Float64

    nlp_local_solver::MathProgBase.AbstractMathProgSolver
    minlp_local_solver::MathProgBase.AbstractMathProgSolver
    mip_solver::MathProgBase.AbstractMathProgSolver

    discretization_var_pick_algo::Any
    discretization_ratio::Any
    discretization_add_partition_method::Any

    presolve_track_time::Bool
    presolve_do_bound_tightening::Bool
    presolve_maxiter::Int
    presolve_tolerance::Float64
    presolve_bound_tightening_method::Any
    presolve_mip_relaxation::Bool
    presolve_mip_timelimit::Float64

    # other options to be added later on
end

function PODSolver(;
    log_level = 1,
    timeout = Inf,
    maxiter = 99,
    rel_gap = 1e-4,
    tolerance = 1e-6,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),

    discretization_var_pick_algo = 0,           # By default pick all variables
    discretization_ratio = 4,
    discretization_add_partition_method = nothing, # Not ready for implementation

    presolve_track_time = false,
    presolve_do_bound_tightening = false,
    presolve_maxiter = 9999,
    presolve_tolerance = 1e-3,
    presolve_bound_tightening_method = 1,
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
    PODSolver(log_level, timeout, maxiter, rel_gap, tolerance,
        deepcopy(nlp_local_solver), deepcopy(minlp_local_solver), deepcopy(mip_solver),
        discretization_var_pick_algo, discretization_ratio, discretization_add_partition_method,
        presolve_track_time, presolve_do_bound_tightening, presolve_maxiter, presolve_tolerance, presolve_bound_tightening_method,presolve_mip_relaxation,presolve_mip_timelimit)
    end

# Create POD nonlinear model: can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)
    if !applicable(MathProgBase.NonlinearModel, s.nlp_local_solver)
        error("NLP local solver $(s.nlp_local_solver) specified is not a NLP solver recognized by POD\n")
    end

    # Translate options into old nonlinearmodel.jl fields
    log_level = s.log_level
    timeout = s.timeout
    maxiter = s.maxiter
    rel_gap = s.rel_gap
    tolerance = s.tolerance
    nlp_local_solver = s.nlp_local_solver
    minlp_local_solver = s.minlp_local_solver
    mip_solver = s.mip_solver
    discretization_var_pick_algo = s.discretization_var_pick_algo
    discretization_ratio = s.discretization_ratio
    discretization_add_partition_method = s.discretization_add_partition_method

    presolve_track_time = s.presolve_track_time
    presolve_do_bound_tightening = s.presolve_do_bound_tightening
    presolve_maxiter = s.presolve_maxiter
    presolve_tolerance = s.presolve_tolerance
    presolve_bound_tightening_method = s.presolve_bound_tightening_method
    presolve_mip_relaxation = s.presolve_mip_relaxation
    presolve_mip_timelimit = s.presolve_mip_timelimit

    return PODNonlinearModel(log_level, timeout, maxiter, rel_gap, tolerance,
                            nlp_local_solver, minlp_local_solver, mip_solver,
                            discretization_var_pick_algo, discretization_ratio, discretization_add_partition_method,
                            presolve_track_time, presolve_do_bound_tightening, presolve_maxiter, presolve_tolerance,presolve_bound_tightening_method,presolve_mip_relaxation,presolve_mip_timelimit)
end

"""
    A small function used to fetch differnt MIP solver setting code
"""
function fetch_timeleft_symbol(m::PODNonlinearModel; kwargs...)
    if string(m.mip_solver)[1:10] == "CPLEX.Cple"
        return :CPX_PARAM_TILIM
    elseif string(m.mip_solver)[1:10] == "Gurobi.Gur"
        return :TimeLimit
    elseif string(m.mip_solver)[1:10] == "Cbc.CbcMat"
        return :seconds
    else found == nothing
        error("Needs support for this MIP solver")
    end
    return
end
