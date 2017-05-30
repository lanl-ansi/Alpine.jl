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

    discrete_vars_choice::Int
    discrete_ratio::Any

    # other options to be added later on
end

function PODSolver(;
    log_level = 1,
    timeout = Inf,
    rel_gap = 1e-5,

    nlp_local_solver = UnsetSolver(),
    minlp_local_solver = UnsetSolver(),
    mip_solver = UnsetSolver(),

    discrete_vars_choice = 0,
    discrete_ratio = 4,
    )

    if nlp_local_solver == UnsetSolver()
        error("No NLP local solver specified (set nlp_local_solver)\n")
    end

    if mip_solver == UnsetSolver()
        error("NO MIP solver specififed (set mip_solver)\n")
    end

    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(log_level, timeout, rel_gap,
        deepcopy(nlp_local_solver), deepcopy(minlp_local_solver), deepcopy(mip_solver),
        discrete_vars_choice, discrete_ratio)
end

# Create POD nonlinear model: can solve with nonlinear algorithm only
function MathProgBase.NonlinearModel(s::PODSolver)
    if !applicable(MathProgBase.NonlinearModel, s.nlp_local_solver)
        error("NLP local solver $(s.nlp_local_solver) specified is not a NLP solver recognized by POD\n")
    end

    # Translate options into old nonlinearmodel.jl fields
    log_level = s.log_level
    timeout = s.timeout
    rel_gap = s.rel_gap
    nlp_local_solver = s.nlp_local_solver
    minlp_local_solver = s.minlp_local_solver
    mip_solver = s.mip_solver
    discrete_vars_choice = s.discrete_vars_choice
    discrete_ratio = s.discrete_ratio

    return PODNonlinearModel(log_level, timeout, rel_gap, nlp_local_solver, minlp_local_solver, mip_solver, discrete_vars_choice, discrete_ratio)
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
end

function set_mip_time_limit(m::PODNonlinearModel)
    for i in 1:length(m.mip_solver.options)
        if fetch_timeleft_symbol(m) in collect(m.mip_solver.options[i])
            deleteat!(m.mip_solver.options, i)
            break
        end
    end
    push!(m.mip_solver.options, (fetch_timeleft_symbol(m), m.timeleft))
end

function print_iis_gurobi(m::Model)

    grb = MathProgBase.getrawsolver(internalmodel(m))
    Gurobi.computeIIS(grb)
    numconstr = Gurobi.num_constrs(grb)
    numvar = Gurobi.num_vars(grb)

    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

    info("Irreducible Inconsistent Subsystem (IIS)")
    info("Variable bounds:")
    for i in 1:numvar
        v = Variable(m, i)
        if iislb[i] != 0 && iisub[i] != 0
            println(getlowerbound(v), " <= ", getname(v), " <= ", getupperbound(v))
        elseif iislb[i] != 0
            println(getname(v), " >= ", getlowerbound(v))
        elseif iisub[i] != 0
            println(getname(v), " <= ", getupperbound(v))
        end
    end

    info("Constraints:")
    for i in 1:numconstr
        if iisconstr[i] != 0
            println(m.linconstr[i])
        end
    end

end
