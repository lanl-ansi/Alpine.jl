using JuMP

type PODNonlinearModel <: MathProgBase.AbstractNonlinearModel
    # solver parameters
    log_level::Int                                              # Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info
    timeout::Float64                                            # Time limit for algorithm (in seconds)
    rel_gap::Float64                                            # Relative optimality gap termination condition
    discrete_vars_choice::Int                                   # [SW] 1: Minimum vertex cover, 0:Max cover

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
    # [SW] is managing this section ----------------------------
    basic_model_mip::JuMP.Model                                 # JuMP convex MIP model for bounding
    iter::Int                                                   # Iteration count
    timeleft::Float64                                           # Parameters used for algorithm
    x_int::Vector{JuMP.Variable}                                # JuMP vector of integer variables (:Int, :Bin)
    x_cont::Vector{JuMP.Variable}                               # JuMP vector of continuous variables
    dict_nonlinear_info::Dict{Any,Any}                          # Dictionary containing details of lifted terms
    lifted_obj_expr_mip::Expr                                   # Lifted objective expression, if linear, same as obj_expr_orig
    lifted_constr_expr_mip::Vector{Expr}                        # Lifted constraints, if linear, same as corresponding constr_expr_orig
    lifted_x_cont::Int                                          # Count of lifted variables
    lifted_obj_aff_mip::Dict{Any, Any}                          # Affine function of the lifted objective expression
    lifted_constr_aff_mip::Vector{Dict{Any, Any}}               # Affine function of the lifted constraints
    discretization::Dict{Any,Any}                               # Discretization points keyed by the variables
    discrete_x_cnt::Int                                         # Number of variables got discretized
    discrete_x::Vector{Any}                                     # Variables on which discretization is performed
    discrete_ratio::Float64                                     # Discretization ratio (use a fixed value for now, later switch to a function)
    sol_incumb_lb::Vector{Float64}                              # Incumbent lower bounding solution
    sol_incumb_ub::Vector{Float64}                              # Incumbent upper bounding solution
    # [SW] is managing this section ----------------------------

    # Solution and bound information
    best_bound::Float64                                         # Best bound from MIP
    best_obj::Float64                                           # Best feasible objective value
    best_sol::Vector{Float64}                                   # Best feasible solution
    gap_rel_opt::Float64                                        # Relative optimality gap = |best_bound - best_obj|/|best_obj|
    final_soln::Vector{Float64}                                 # Final solution

    # Logging information and status
    logs::Dict{Symbol,Any}                                      # Logging information
    status::Symbol                                              # Current POD status
    lb_status::Symbol                                           # [SW] Lower bound status

    # constructor
    function PODNonlinearModel(log_level, timeout, rel_gap, nlp_local_solver, minlp_local_solver, mip_solver)
        m = new()
        m.log_level = log_level
        m.timeout = timeout
        m.rel_gap = rel_gap
        m.discrete_vars_choice = 0          #[SW] Added default value on how to choose discretize variables

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
        m.dict_nonlinear_info = Dict()

        m.lifted_constr_expr_mip = []       # [SW] added
        m.lifted_constr_aff_mip = []        # [SW] added
        m.discrete_x = []                   # [SW] added
        m.discretization = Dict()           # [SW] added
        m.discrete_ratio = 4                # [SW] added
        m.timeleft = m.timeout              # [SW] added
        m.iter = 0                          # [SW] added

        m.best_obj = Inf
        m.best_bound = -Inf
        m.gap_rel_opt = NaN

        m.status = :NotLoaded
        m.lb_status = :NotLoaded
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
    # m.var_type_orig = fill(:Cont, m.num_var_orig)
    m.var_type_orig = [getcategory(Variable(d.m, i)) for i in 1:m.num_var_orig] #[SW] getting ready for MINLP

    # Summarize constraints information in original model
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

    # [SW] POD_MIP model
    populate_dict_nonlinear_info(m)
    populate_lifted_expr(m)
    m.lifted_x_cont = length(m.dict_nonlinear_info)
    populate_lifted_affine(m)

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
    global_solve(m)
end

MathProgBase.setwarmstart!(m::PODNonlinearModel, x) = (m.var_start_orig = x)
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
    initialize_discretization(m) # [SW] initialize discretization here
end

#=
    Main Adaptive Partitioning (AP) Algorithm.
=#
function global_solve(m::PODNonlinearModel;kwargs...)
    while (m.best_obj - m.best_bound)/m.best_obj < m.gap_rel_opt || m.timeleft > 0.0
        m.iter += 1
        set_mip_time_limit(m)       # Regulates timing
        lower_bounding_mip(m)       # Build lower bounding model
        lb_solve(m)                 # Solve lower bounding model
        # print(m.basic_model_mip)
        add_discretization(m)       # Manage discretization
        ub_solve(m)                 # Solve upper bounding model
        m.rel_gap = (m.best_obj - m.best_bound)/m.best_obj
        @printf "gap updates %13.5f.\n" m.rel_gap
        if m.iter > 5
            error("TESTING...")
        end
    end
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
    if local_solve_nlp_status == :Optimal || local_solve_nlp_status == :Suboptimal || local_solve_nlp_status == :UserLimit  # [SW] Added :UserLimit for feasibility tuning
        m.best_sol = MathProgBase.getsolution(local_solve_nlp_model)
        m.sol_incumb_ub = copy(m.best_sol)  # [SW] added for temp holder
        m.best_obj = MathProgBase.getobjval(local_solve_nlp_model)
        m.status = local_solve_nlp_status
        (m.log_level > 0) && @printf "best feasible solution objective = %13.5f.\n" m.best_obj
        return
    elseif local_solve_nlp_status == :Infeasible
        (m.log_level > 0) && println("problem infeasible.")
        m.status = :Infeasible
        return
        # TODO Figure out the conditions for this to hold!
        # [SW] Past experience, given different solver utilized,
        # BOMIN -> Resolve multiple times at root relaxation will help
        # KNITRO -> Use multi-start will help
        # Ipopt -> Unknown method
        # This may need to be conducted given specific NLP solver utilized
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

function ub_solve(m::PODNonlinearModel; kwargs...)

    convertor = Dict(:Max=>:>, :Min=>:<)

    ub_solve_nlp_model = MathProgBase.NonlinearModel(m.nlp_local_solver)
    l_var, u_var = tight_ub_bounds(m)
    MathProgBase.loadproblem!(ub_solve_nlp_model, m.num_var_orig, m.num_constr_orig, l_var, u_var,m.l_constr_orig, m.u_constr_orig, m.sense_orig, m.d_orig)
    MathProgBase.setwarmstart!(ub_solve_nlp_model, m.best_sol[1:m.num_var_orig])

    start_ub_solve = time()
    MathProgBase.optimize!(ub_solve_nlp_model)
    cputime_ub_solve = time() - start_ub_solve
    m.logs[:total_time] += cputime_ub_solve
    m.timeleft -= cputime_ub_solve

    status = MathProgBase.status(ub_solve_nlp_model)
    if status == :Optimal || status == :Suboptimal || status == :UserLimit  # [SW] Added :UserLimit for feasibility tuning
        candidate_obj = MathProgBase.getobjval(ub_solve_nlp_model)
        if eval(convertor[m.sense_orig])(candidate_obj, m.best_obj)
            m.best_sol = MathProgBase.getsolution(ub_solve_nlp_model)
            m.sol_incumb_ub = copy(m.best_sol)  # [SW] added for temp holder
            (m.log_level > 0) && @printf "[+] best feasible solution objective = %13.5f.\n" m.best_obj
        else
            @printf "[-] best feasible solution objective = %13.5f.\n" m.best_obj
        end
    elseif ub_solve_nlp_status == :Infeasible
        (m.log_level > 0) && println("Incumbent unchanged due to infeasible UB solve.")
    elseif ub_solve_nlp_status == :Unbounded
        error("ERROR happend. Unbounded problem during UB solve")
    else
        error("NLP solver failure during a UB solve.")
    end
end

function lb_solve(m::PODNonlinearModel; kwargs...)

    convertor = Dict(:Max=>:<, :Min=>:>)
    set_mip_time_limit(m)
    start_lb_solve = time()
    status = solve(m.basic_model_mip, suppress_warnings=true)
    cputime_lb_solve = time() - start_lb_solve
    m.logs[:total_time] += cputime_lb_solve
    m.timeleft -= cputime_lb_solve
    if status == :Optimal || status == :Suboptimal || status == :UserLimit
        candidate_bound = getobjectivevalue(m.basic_model_mip)
        if eval(convertor[m.sense_orig])(candidate_bound, m.best_bound)
            m.best_bound = candidate_bound
            m.sol_incumb_lb = [getvalue(Variable(m.basic_model_mip, i)) for i in 1:m.num_var_orig]
            m.lb_status = status
            (m.log_level > 0) && @printf "[+] LB incumbent = %13.5f.\n" m.best_bound
        end
    elseif status == :Infeasible
        print_iis_gurobi(m.basic_model_mip)
        #= Case when this happens::
            1) :UserLimits ? but will this return :UserLimits?
            2) :Numerical Difficult
        =#
        error("[MIP INFEASIBLE] There is some issue about LB problem")
    elseif status == :Unbounded
        error("[MIP UNBOUNDED] MIP solver failure")
    else
        error("[MIP UNKNOWN] MIP solver failure.")
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

# Create dictionary of statuses for POD algorithm
function create_status!(m)
    status = Dict{Symbol,Symbol}()

    status[:local_solve] = :none            # Status of local solve
    status[:bounding_solve] = :none         # Status of bounding solve
    status[:bnd_tightening_solve] = :none   # Status of bound-tightening solve
    status[:pod] = :none                    # Status of POD algorithm

    m.status = status
end
