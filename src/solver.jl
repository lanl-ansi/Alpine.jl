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

    var_discretization_algo::Int
    discretization_ratio::Any
    method_pick_vars_discretization::Function

    # other options to be added later on
end

function PODSolver(;
    log_level = 1,
    timeout = Inf,
    maxiter = 99,
    rel_gap = 1e-4,
    tolerance = 1e-4,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),

    var_discretization_algo = 0,
    discretization_ratio = 4,
    method_pick_vars_discretization = min_vertex_cover,
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
        var_discretization_algo, discretization_ratio,
        method_pick_vars_discretization)
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
    var_discretization_algo = s.var_discretization_algo
    discretization_ratio = s.discretization_ratio
    method_pick_vars_discretization = s.method_pick_vars_discretization

    return PODNonlinearModel(log_level, timeout, maxiter, rel_gap, tolerance,
                            nlp_local_solver, minlp_local_solver, mip_solver,
                            var_discretization_algo, discretization_ratio,
                            method_pick_vars_discretization)
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

function update_time_limit(m::PODNonlinearModel)
    for i in 1:length(m.mip_solver.options)
        if fetch_timeleft_symbol(m) in collect(m.mip_solver.options[i])
            deleteat!(m.mip_solver.options, i)
            break
        end
    end
    if m.timeout != Inf
        push!(m.mip_solver.options, (fetch_timeleft_symbol(m), max(0.0,m.timeout - m.logs[:total_time])))
    end
    return
end
