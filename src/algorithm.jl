using JuMP

type PODNonlinearModel <: MathProgBase.AbstractNonlinearModel
    # solver parameters
    verbose::Int
    time_limit::Float64
    opt_tolerance::Float64

    # add all the solver options
    nlp_local_solver::MathProgBase.AbstractMathProgSolver
    minlp_local_solver::MathProgBase.AbstractMathProgSolver
    mip_solver::MathProgBase.AbstractMathProgSolver

    # other options go here

    # initial data
    num_var::Int
    num_int_var::Int 
    num_constr::Int
    num_nlconstr::Int
    solution::Vector{Float64}
    status
    objval::Float64
    objbound::Float64
    iterations::Int

    A
    A_lb
    A_ub
    lb::Vector{Float64}
    ub::Vector{Float64}
    l::Vector{Float64}
    u::Vector{Float64}
    c::Vector{Float64}
    vartype::Vector{Symbol}
    constrtype::Vector{Symbol}
    constrlinear::Vector{Bool}
    objsense::Symbol
    objlinear::Bool
    d # jacobian for NL solve

    nlp_load_timer
    total_time

    # constructor
    function PODNonlinearModel(verbose, time_limit, opt_tolerance, nlp_local_solver, minlp_local_solver, mip_solver)
        m = new()
        m.verbose = verbose
        m.time_limit = time_limit
        m.opt_tolerance = opt_tolerance
        m.nlp_local_solver = nlp_local_solver
        m.minlp_local_solver = minlp_local_solver
        m.mip_solver = mip_solver
        m.total_time = 0.
        return m
    end
end

function MathProgBase.loadproblem!(m::PODNonlinearModel, numVar, numConstr, l, u, lb, ub, sense, d)

    #if !applicable(MathProgBase.NonlinearModel, m.nlp_local_solver)
    #    error("$(m.nlp_local_solver) is not a nonlinear solver.")
    #end
    m.num_var = numVar
    m.num_constr = numConstr
    m.lb = lb
    m.ub = ub
    m.l = l
    m.u = u
    m.objsense = sense
    m.d = d
    m.vartype = fill(:Cont,numVar)

    MathProgBase.initialize(d, [:Grad,:Jac,:Hess])

    m.constrtype = Array(Symbol, numConstr)
    for i = 1:numConstr
        if lb[i] > -Inf && ub[i] < Inf
            m.constrtype[i] = :(==)
        elseif lb[i] > -Inf
            m.constrtype[i] = :(>=)
        else
            m.constrtype[i] = :(<=)
        end
    end
    m.solution = fill(NaN, m.num_var)
end

function populatelinearmatrix(m::PODNonlinearModel)
    # set up map of linear rows
    constrlinear = Array(Bool, m.num_constr)
    numlinear = 0
    constraint_to_linear = fill(-1,m.num_constr)
    for i = 1:m.num_constr
        constrlinear[i] = MathProgBase.isconstrlinear(m.d, i)
        if constrlinear[i]
            numlinear += 1
            constraint_to_linear[i] = numlinear
        end
    end
    m.num_nlconstr = m.num_constr - numlinear

    # extract sparse jacobian structure
    jac_I, jac_J = MathProgBase.jac_structure(m.d)

    # evaluate jacobian at x = 0
    c = zeros(m.num_var)
    x = m.solution
    jac_V = zeros(length(jac_I))
    MathProgBase.eval_jac_g(m.d, jac_V, x)
    MathProgBase.eval_grad_f(m.d, c, x)
    m.objlinear = MathProgBase.isobjlinear(m.d)
    if m.objlinear
        (m.verbose > 0) && println("Objective function is linear")
        m.c = c
    else
        (m.verbose > 0) && println("Objective function is nonlinear")
        m.c = zeros(m.num_var)
    end

    # Build up sparse matrix for linear constraints
    A_I = Int[]
    A_J = Int[]
    A_V = Float64[]

    for k in 1:length(jac_I)
        row = jac_I[k]
        if !constrlinear[row]
            continue
        end
        row = constraint_to_linear[row]
        push!(A_I,row); push!(A_J, jac_J[k]); push!(A_V, jac_V[k])
    end

    m.A = sparse(A_I, A_J, A_V, m.num_constr-m.num_nlconstr, m.num_var)

    # g(x) might have a constant, i.e., a'x + b
    # let's find b
    constraint_value = zeros(m.num_constr)
    MathProgBase.eval_g(m.d, constraint_value, x)
    b = constraint_value[constrlinear] - m.A * x
    # so linear constraints are of the form lb ≤ a'x + b ≤ ub

    # set up A_lb and A_ub vectors
    m.A_lb = m.lb[constrlinear] - b
    m.A_ub = m.ub[constrlinear] - b

    # Now we have linear parts
    m.constrlinear = constrlinear
    return
end

function MathProgBase.optimize!(m::PODNonlinearModel)
    start = time()

    cputime_nlp = 0.0
    cputime_mip = 0.0

    m.nlp_load_timer = 0.0

    # if we haven't gotten a starting point,
    # (e.g., if acting as MIQP solver) assume zero is okay
    if any(isnan,m.solution)
        m.solution = zeros(length(m.solution))
    end

    populatelinearmatrix(m)
    #MathProgBase.initialize(m.d, [:Grad, :Jac, :ExprGraph, :Hess, :HessVec])
    jac_I, jac_J = MathProgBase.jac_structure(m.d)
    jac_V = zeros(length(jac_I))
    grad_f = zeros(m.num_var)
    
    ini_nlp_model = MathProgBase.NonlinearModel(m.nlp_local_solver)
    start_load = time()
    MathProgBase.loadproblem!(ini_nlp_model, m.num_var, m.num_constr, m.l, m.u, m.lb, m.ub, m.objsense, m.d)
    m.nlp_load_timer += time() - start_load

    # pass in starting point
    #MathProgBase.setwarmstart!(nlp_model, m.solution)
    MathProgBase.setwarmstart!(ini_nlp_model, m.solution[1:m.num_var])

    start_nlp = time()
    MathProgBase.optimize!(ini_nlp_model)
    cputime_nlp += time() - start_nlp

    ini_nlp_status = MathProgBase.status(ini_nlp_model)
    if ini_nlp_status == :Optimal || ini_nlp_status == :Suboptimal
        m.solution = MathProgBase.getsolution(ini_nlp_model)
        m.objval = MathProgBase.getobjval(ini_nlp_model)
        m.status = ini_nlp_status
        return
    elseif ini_nlp_status == :Infeasible
        (m.verbose > 0) && println("Initial NLP Infeasible.")
        m.status = :Infeasible
        return
        # TODO Figure out the conditions for this to hold!
    elseif ini_nlp_status == :Unbounded
        warn("Initial NLP Relaxation Unbounded.")
        m.status = :Unbounded
        return
    else
        warn("NLP Solver Failure.")
        m.status = :Error
        return
    end
    ini_nlp_objval = MathProgBase.getobjval(ini_nlp_model)

    (m.verbose > 0) && println("\nPOD started...\n")
    (m.verbose > 0) && println("MINLP algorithm $(m.algorithm) is chosen.")
    (m.verbose > 0) && println("MINLP has $(m.num_var) variables, $(m.num_constr - m.num_nlconstr) linear constraints, $(m.num_nlconstr) nonlinear constraints.")
    (m.verbose > 0) && @printf "Initial objective = %13.5f.\n\n" ini_nlp_objval

    m.status = :UserLimit
    m.objval = Inf

    m.total_time = time()-start

    if m.verbose > 0
        println("\nPOD finished...\n")
        @printf "Status            = %13s.\n" m.status
        (m.status == :Optimal) && @printf "Optimum objective = %13.5f.\n" m.objval
        @printf "Total time        = %13.5f sec.\n" m.total_time
        @printf "NLP total time    = %13.5f sec.\n" cputime_nlp
        @printf "Subprob load time = %13.5f sec.\n" m.nlp_load_timer
    end
end

MathProgBase.setwarmstart!(m::PODNonlinearModel, x) = (m.solution = x)
MathProgBase.setvartype!(m::PODNonlinearModel, v::Vector{Symbol}) = (m.vartype = v)
#MathProgBase.setvarUB!(m::IpoptMathProgModel, v::Vector{Float64}) = (m.u = v)
#MathProgBase.setvarLB!(m::IpoptMathProgModel, v::Vector{Float64}) = (m.l = v)

MathProgBase.status(m::PODNonlinearModel) = m.status
MathProgBase.getobjval(m::PODNonlinearModel) = m.objval
MathProgBase.getobjbound(m::PODNonlinearModel) = m.objbound
MathProgBase.getsolution(m::PODNonlinearModel) = m.solution
MathProgBase.getsolvetime(m::PODNonlinearModel) = m.total_time
