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
    constr_expr_orig::Vector{Expr}                              # Constraint expressions
    obj_expr_orig::Expr                                         # Objective expression

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
    num_lconstr_updated::Int                                    # Updated number of linear constraints - includes linear constraints added via @NLconstraint macro
    num_nlconstr_updated::Int                                   # Updated number of non-linear constraints
    indexes_lconstr_updated::Vector{Int}                        # Indexes of updated linear constraints

    # local solution model extra data for each iteration
    l_var::Vector{Float64}                                      # Updated variable lower bounds for local solve
    u_var::Vector{Float64}                                      # Updated variable upper bounds for local solve
    var_type::Vector{Symbol}                                    # Updated variable type for local solve

    # mixed-integer convex program bounding model
    model_mip::JuMP.Model                                       # JuMP convex MIP model for bounding
    x_int::Vector{JuMP.Variable}                                # JuMP vector of integer variables (:Int, :Bin)
    x_cont::Vector{JuMP.Variable}                               # JuMP vector of continuous variables
    dict_nonlinear_terms::Dict{Any,Any}                          # Dictionary containing details of lifted terms
    var_discretization::Vector{Any}                             # Variables on which discretization is performed
    discretization::Dict{Any,Set{Float64}}                      # Discretization points keyed by the variables
    lifted_obj_expr_mip::Expr                                   # Lifted objective expression, if linear, same as obj_expr_orig [SW]
    lifted_constr_expr_mip::Vector{Expr}                         # Lifted constraints, if linear, same as corresponding constr_expr_orig [SW]


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
        m.constr_expr_orig = Expr[]
        m.num_lconstr_updated = 0
        m.num_nlconstr_updated = 0
        m.indexes_lconstr_updated = Int[]
        # Need to be initialized [SW]
        m.dict_nonlinear_terms = Dict()


        m.best_obj = Inf
        m.best_bound = -Inf
        m.gap_rel_opt = NaN

        m.status = :NotLoaded
        #create_status!(m)
        create_logs!(m)


        return m
    end
end

function MathProgBase.loadproblem!(m::PODNonlinearModel,
                                    num_var::Int, num_constr::Int,
                                    l_var::Vector{Float64}, u_var::Vector{Float64},
                                    l_constr::Vector{Float64}, u_constr::Vector{Float64},
                                    sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

    m.num_var_orig = num_var
    m.num_constr_orig = num_constr
    m.l_var_orig = l_var
    m.u_var_orig = u_var
    m.l_constr_orig = l_constr
    m.u_constr_orig = u_constr
    m.sense_orig = sense
    if m.sense_orig == :Max
        m.best_obj = -Inf
        m.best_bound = Inf
    end
    m.d_orig = d
    MathProgBase.initialize(m.d_orig, [:Grad,:Jac,:Hess,:ExprGraph])
    for i in 1:m.num_constr_orig
        push!(m.constr_expr_orig, MathProgBase.constr_expr(d, i))
    end
    m.obj_expr_orig = MathProgBase.obj_expr(d)

    

    # this should change later for an MINLP
    m.var_type_orig = fill(:Cont, m.num_var_orig)

    m.constr_type_orig = Array(Symbol, m.num_constr_orig)
    for i in 1:m.num_constr_orig
        if l_constr[i] > -Inf && u_constr[i] < Inf
            m.constr_type_orig[i] = :(==)
        elseif l_constr[i] > -Inf
            m.constr_type_orig[i] = :(>=)
        else
            m.constr_type_orig[i] = :(<=)
        end
    end

    for i = 1:m.num_constr_orig
        if MathProgBase.isconstrlinear(m.d_orig, i)
            m.num_lconstr_orig += 1
        end
    end
    m.num_nlconstr_orig = m.num_constr_orig - m.num_lconstr_orig


    # not using this any where (in optional fields)
    m.is_obj_linear_orig = MathProgBase.isobjlinear(m.d_orig)

    m.best_sol = fill(NaN, m.num_var_orig)

    if m.log_level > 0
        @printf "full problem loaded into POD.\n"
        @printf "number of constraints = %d.\n" m.num_constr_orig
        @printf "number of non-linear constraints = %d.\n" m.num_nlconstr_orig
        @printf "number of linear constraints = %d.\n" m.num_lconstr_orig
        @printf "number of variables = %d.\n" m.num_var_orig
    end
end

function MathProgBase.optimize!(m::PODNonlinearModel)
    start = time()
    if any(isnan, m.best_sol)
        m.best_sol = zeros(length(m.best_sol))
    end

    presolve(m)
    # presolve(m)
    # bounding tightening
    # bounding solve
    # compute gap and update the gap

    #=
    while (gap > tol)
        local_solve(m)
        bounding_solve(m)
        update_gap(m)
    end
    =#

    # create a bounding convex mip model (JuMP)
end

MathProgBase.setwarmstart!(m::PODNonlinearModel, x) = (m.best_sol = x)
MathProgBase.setvartype!(m::PODNonlinearModel, v::Vector{Symbol}) = (m.var_type_orig = v)

MathProgBase.status(m::PODNonlinearModel) = m.status
MathProgBase.getobjval(m::PODNonlinearModel) = m.best_obj
MathProgBase.getobjbound(m::PODNonlinearModel) = m.best_bound
MathProgBase.getsolution(m::PODNonlinearModel) = m.best_sol
MathProgBase.getsolvetime(m::PODNonlinearModel) = m.logs[:total_time]

function presolve(m::PODNonlinearModel)
    (m.log_level > 0) && println("\nPOD algorithm presolver started.")
    (m.log_level > 0) && println("1. performing local solve to obtain a feasible solution.")
    local_solve(m)
end

function local_solve(m::PODNonlinearModel)
    local_solve_nlp_model = MathProgBase.NonlinearModel(m.nlp_local_solver)
    MathProgBase.loadproblem!(local_solve_nlp_model, m.num_var_orig, m.num_constr_orig, m.l_var_orig, m.u_var_orig, m.l_constr_orig, m.u_constr_orig, m.sense_orig, m.d_orig)
    MathProgBase.setwarmstart!(local_solve_nlp_model, m.best_sol[1:m.num_var_orig])

    start_local_solve = time()
    MathProgBase.optimize!(local_solve_nlp_model)
    cputime_local_solve = time() - start_local_solve
    m.logs[:total_time] += cputime_local_solve

    local_solve_nlp_status = MathProgBase.status(local_solve_nlp_model)
    if local_solve_nlp_status == :Optimal || local_solve_nlp_status == :Suboptimal
        m.best_sol = MathProgBase.getsolution(local_solve_nlp_model)
        m.best_obj = MathProgBase.getobjval(local_solve_nlp_model)
        m.status = local_solve_nlp_status
        (m.log_level > 0) && @printf "best feasible solution objective = %13.5f.\n" m.best_obj
        return
    elseif local_solve_nlp_status == :Infeasible
        (m.log_level > 0) && println("problem infeasible.")
        m.status = :Infeasible
        return
        # TODO Figure out the conditions for this to hold!
    elseif local_solve_nlp_status == :Unbounded
        warn("initial problem unbounded.")
        m.status = :Unbounded
        return
    else
        warn("NLP solver failure.")
        m.status = :Error
        return
    end
    local_solve_nlp_objval = MathProgBase.getobjval(local_solve_nlp_model)

    m.best_obj = local_solve_nlp_objval

    return

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

# Create dictionary of statuses for POD algorithm
function create_status!(m)
    status = Dict{Symbol,Symbol}()

    status[:local_solve] = :none            # Status of local solve
    status[:bounding_solve] = :none         # Status of bounding solve
    status[:bnd_tightening_solve] = :none   # Status of bound-tightening solve
    status[:pod] = :none                    # Status of POD algorithm

    m.status = status
end
