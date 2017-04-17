export PODSolver

# Dummy solver
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

type PODSolver <: MathProgBase.AbstractMathProgSolver 
    log_level::Int
    timeout::Float64
    rel_gap::Float64

    mip_solver::MathProgBase.AbstractMathProgSolver
    cont_solver::MathProgBase.AbstractMathProgSolver 
end 

function PODSolver(;
    log_level = 1,
    timeout = Inf,
    rel_gap = 1e-5,

    mip_solver = UnsetSolver(),
    cont_solver = UnsetSolver())


    # Deepcopy the solvers because we may change option values inside POD
    PODSolver(log_level, timeout, rel_gap, deepcopy(mip_solver), deepcopy(cont_solver))
end


