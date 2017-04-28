using JuMP

type PODNonlinearModel <: MathProgBase.AbstractNonlinearModel
    # solver parameters
    log_level::Int                                              # Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info
    timeout::Float64                                            # Time limit for algorithm (in seconds)
    rel_gap::Float64                                            # Relative optimality gap termination condition

    # add all the solver options
    nlp_local_solver::MathProgBase.AbstractMathProgSolver       # Local continuous NLP solver for solving NLPs at each iteration
    minlp_local_solver::MathProgBase.AbstractMathProgSolver     # Local MINLP solver for solving MINLPs at each iteration
    mip_solver::MathProgBase.AbstractMathProgSolver             # MILP solver for successive lower bound solves

    # other options go here

    # initial data provided by user
    num_var_orig::Int                                           # Initial number of variables
    num_cont_var_orig::Int                                      # Initial number of continuous variables
    num_int_var_orig::Int                                       # Initial number of binary/integer variables
    num_constr_orig::Int                                        # Initial number of constraints
    num_lconstr_orig::Int                                       # Initial number of linear constraints
    num_nlconstr_orig::Int                                      # Initial number of non-linear constraints
    var_type_orig::Vector{Symbol}                               # Variable type vector on original variables (only :Bin, :Cont, :Int)
    var_start_orig::Vector{Float64}                             # Variable warm start vector on original variables
    constr_type_orig::Vector{Symbol}                            # Constraint type vector on original variables (only :(==), :(>=), :(<=))

    # extra initial data that is useful for non-linear local continuous solves
    l_var_orig::Vector{Float64}                                 # Variable lower bounds
    u_var_orig::Vector{Float64}                                 # Variable upper bounds
    l_constr_orig::Vector{Float64}                              # Constraint lower bounds
    u_constr_orig::Vector{Float64}                              # Constraint upper bounds
    sense_orig::Symbol                                          # Problem type (:Min, :Max)
    d_orig::JuMP.NLPEvaluator                                   # Instance of AbstractNLPEvaluator for evaluating gradient, Hessian-vector products, and Hessians of the Lagrangian

    # additional initial data that may be useful later on - non populated for now
    A_orig                                                      # Linear constraint matrix
    A_l_orig                                                    # Linear constraint matrix LHS
    A_u_orig                                                    # Linear constraint matrix RHS
    is_obj_linear_orig::Bool                                    # Bool variable for type of objective
    c_orig::Vector{Float64}                                     # Coefficient vector for linear objective

    # local solution model extra data for each iteration
    l_var::Vector{Float64}                                      # Updated variable lower bounds for local solve
    u_var::Vector{Float64}                                      # Updated variable upper bounds for local solve

    # mixed-integer convex program bounding model
    model_mip::JuMP.Model                                       # JuMP convex MIP model for bounding
    x_int::Vector{JuMP.Variable}                                # JuMP vector of integer variables (:Int, :Bin)
    x_cont::Vector{JuMP.Variable}                               # JuMP vector of continuous variables

    # Solution and bound information
    best_bound::Float64                                         # Best bound from MIP
    best_obj::Float64                                           # Best feasible objective value
    best_sol::Vector{Float64}                                   # Best feasible solution
    gap_rel_opt::Float64                                        # Relative optimality gap = |best_bound - best_obj|/|best_obj|
    final_soln::Vector{Float64}                                 # Final solution

    # Logging information and status
    logs::Dict{Symbol,Any}                                      # Logging information
    status::Symbol                                              # Current POD status

    # constructor
    function PODNonlinearModel(log_level, timeout, rel_gap, nlp_local_solver, minlp_local_solver, mip_solver)
        m = new()
        m.log_level = log_level
        m.timeout = timeout
        m.rel_gap = rel_gap

        m.nlp_local_solver = nlp_local_solver
        m.minlp_local_solver = minlp_local_solver
        m.mip_solver = mip_solver

        m.num_var_orig = 0
        m.num_cont_var_orig = 0
        m.num_int_var_orig = 0
        m.num_constr_orig = 0
        m.num_lconstr_orig = 0
        m.num_nlconstr_orig = 0
        m.var_type_orig = Symbol[]
        m.var_start_orig = Float64[]
        m.constr_type_orig = Symbol[]

        m.best_obj = Inf
        m.best_bound = -Inf
        m.gap_rel_opt = NaN
        m.status = :NotLoaded

        create_logs!(m)

        return m
    end
end

#=========================================================
 Logging and printing functions
=========================================================#

# Create dictionary of logs for timing and iteration counts
function create_logs!(m)
    logs = Dict{Symbol,Any}()

    # Timers
    logs[:total_time] = 0.  # Total run-time of the algorithm

    # Counters
    logs[:n_iter] = 0       # Number of iterations in iterative
    logs[:n_feas] = 0       # Number of times get a new feasible solution

    m.logs = logs
end

function MathProgBase.loadproblem!(m:PODNonlinearModel, numVar, numConstr, l, u, lb, ub, sense, d)

    if !applicable(MathProgBase.NonlinearModel, m.nlp_local_solver)
        error("$(m.nlp_local_solver) is not a nonlinear solver.")
    end

    m.num_var_orig = numVar
    m.num_constr_orig = numConstr
    m.l_var_orig = l
    m.u_var_orig = u
    m.l_var = l
    m.u_var = u
    m.l_constr_orig = lb
    m.u_constr_orig = ub
    m.sense_orig = sense
    m.d_orig = d
    # this should change later for an MINLP
    m.var_type_orig = fill(:Cont, m.num_var_orig)

    m.constr_type_orig = Array(Symbol, m.num_constr_orig)
    for i = 1:m.num_constr_orig
        if lb[i] > -Inf && ub[i] < Inf
            m.constr_type_orig[i] = :(==)
        elseif lb[i] > -Inf
            m.constr_type_orig[i] = :(>=)
        else
            m.constr_type_orig[i] = :(<=)
        end
    end

    # not using this any where (in optional fields)
    m.is_obj_linear_orig = MathProgBase.isobjlinear(m.d_orig)

    m.best_sol = fill(NaN, m.num_var_orig)
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
        (m.log_level > 0) && println("Objective function is linear")
        m.c = c
    else
        (m.log_level > 0) && println("Objective function is nonlinear")
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
        (m.log_level > 0) && println("Initial NLP Infeasible.")
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

    (m.log_level > 0) && println("\nPOD started...\n")
    (m.log_level > 0) && println("MINLP algorithm $(m.algorithm) is chosen.")
    (m.log_level > 0) && println("MINLP has $(m.num_var) variables, $(m.num_constr - m.num_nlconstr) linear constraints, $(m.num_nlconstr) nonlinear constraints.")
    (m.log_level > 0) && @printf "Initial objective = %13.5f.\n\n" ini_nlp_objval

    m.status = :UserLimit
    m.objval = Inf

    m.total_time = time()-start

    if m.log_level > 0
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
