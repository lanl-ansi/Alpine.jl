export PODSolver

# Dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

type PODSolver <: MathProgBase.AbstractMathProgSolver
    log_level::Int
    timeout::Float64
    rel_gap::Float64

    nlp_local_solver::MathProgBase.AbstractMathProgSolver
    minlp_local_solver::MathProgBase.AbstractMathProgSolver
    mip_solver::MathProgBase.AbstractMathProgSolver
end

function PODSolver(;
    log_level = 1,
    timeout = Inf,
    rel_gap = 1e-5,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),
    )

    if nlp_local_solver == UnsetSolver()
        error("No NLP local solver specified (set nlp_local_solver)\n")
    end

    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(log_level, timeout, rel_gap, deepcopy(nlp_local_solver), deepcopy(minlp_local_solver), deepcopy(mip_solver))
end

# Create POD nonlinear model: can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)
    if !applicable(MathProgBase.NonlinearModel, s.nlp_local_solver)
        error("NLP local solver (nlp_local_solver) specified is not a NLP solver recognized by MathProgBase\n")
    end

    # Translate options into old nonlinearmodel.jl fields
    verbose = s.log_level
    time_limit = s.timeout
    opt_tolerance = s.rel_gap
    nlp_local_solver = s.nlp_local_solver
    minlp_local_solver = s.minlp_local_solver
    mip_solver = s.mip_solver

    return PODNonlinearModel(verbose, time_limit, opt_tolerance, nlp_local_solver, minlp_local_solver, mip_solver)
end
