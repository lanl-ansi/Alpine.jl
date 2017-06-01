using JuMP

type PODNonlinearModel <: MathProgBase.AbstractNonlinearModel
    # solver parameters
    log_level::Int                                              # Verbosity flag: 0 for quiet, 1 for basic solve info, 2 for iteration info
    timeout::Float64                                            # Time limit for algorithm (in seconds)
    rel_gap::Float64                                            # Relative optimality gap termination condition
    var_discretization_algo::Int                                # Algorithm for choosing the variables to discretize: 1 for minimum vertex cover, 0 for all variables
    discretization_ratio::Float64                               # Discretization ratio parameter (use a fixed value for now, later switch to a function)

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
    nonlinear_info::Dict{Any,Any}                               # Dictionary containing details of lifted terms
    lifted_obj_expr_mip::Expr                                   # Lifted objective expression; if linear, same as obj_expr_orig
    lifted_constr_expr_mip::Vector{Expr}                        # Lifted constraints; if linear, same as corresponding constr_expr_orig
    num_var_lifted_mip::Int                                     # Number of lifted variables
    lifted_obj_aff_mip::Dict{Any, Any}                          # Lifted objective expression in affine form
    lifted_constr_aff_mip::Vector{Dict{Any, Any}}               # Lifted constraint expressions in affine form
    discretization::Dict{Any,Any}                               # Discretization points keyed by the variables
    num_var_discretization_mip::Int                             # Number of variables on which discretization is performed
    var_discretization_mip::Vector{Any}                         # Variables on which discretization is performed
    sol_incumb_lb::Vector{Float64}                              # Incumbent lower bounding solution
    sol_incumb_ub::Vector{Float64}                              # Incumbent upper bounding solution

    # Solution and bound information
    best_bound::Float64                                         # Best bound from MIP
    best_obj::Float64                                           # Best feasible objective value
    best_sol::Vector{Float64}                                   # Best feasible solution
    best_bound_sol::Vector{Float64}                             # Best bound solution
    best_rel_gap::Float64                                       # Relative optimality gap = |best_bound - best_obj|/|best_obj|
    final_soln::Vector{Float64}                                 # Final solution

    # Logging information and status
    logs::Dict{Symbol,Any}                                      # Logging information
    status::Dict{Symbol,Symbol}                                 # Detailed status of each different phases in algorithm
    pod_status::Symbol                                          # Current POD status

    # constructor
    function PODNonlinearModel(log_level, timeout, rel_gap, nlp_local_solver, minlp_local_solver, mip_solver, var_discretization_algo, discretization_ratio)
        m = new()
        m.log_level = log_level
        m.timeout = timeout
        m.rel_gap = rel_gap
        m.var_discretization_algo = var_discretization_algo

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

        m.nonlinear_info = Dict()
        m.lifted_constr_expr_mip = []
        m.lifted_constr_aff_mip = []
        m.var_discretization_mip = []
        m.discretization = Dict()
        m.discretization_ratio = discretization_ratio

        m.best_obj = Inf
        m.best_bound = -Inf
        m.best_rel_gap = Inf
        m.pod_status = :NotLoaded

        create_status!(m)
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

    m.var_type_orig = [getcategory(Variable(d.m, i)) for i in 1:m.num_var_orig]

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

    # populate data to create the bounding model
    populate_nonlinear_info(m)
    populate_lifted_expr(m)
    m.num_var_lifted_mip = length(m.nonlinear_info)
    populate_lifted_affine(m)

    m.best_sol = fill(NaN, m.num_var_orig)

    logging_summary(m)
end

function MathProgBase.optimize!(m::PODNonlinearModel)
    start = time()
    if any(isnan, m.best_sol)
        m.best_sol = zeros(length(m.best_sol))
    end
    presolve(m)
    global_solve(m)
    summary_status(m)
end

MathProgBase.setwarmstart!(m::PODNonlinearModel, x) = (m.var_start_orig = x)
MathProgBase.setvartype!(m::PODNonlinearModel, v::Vector{Symbol}) = (m.var_type_orig = v)

MathProgBase.status(m::PODNonlinearModel) = m.pod_status
MathProgBase.getobjval(m::PODNonlinearModel) = m.best_obj
MathProgBase.getobjbound(m::PODNonlinearModel) = m.best_bound
MathProgBase.getsolution(m::PODNonlinearModel) = m.best_sol
MathProgBase.getsolvetime(m::PODNonlinearModel) = m.logs[:total_time]

function presolve(m::PODNonlinearModel)

    start_presolve = time()
    (m.log_level > 0) && println("\nPOD algorithm presolver started.")
    (m.log_level > 0) && println("Performing local solve to obtain a feasible solution.")
    local_solve(m, presolve = true)

    if m.status[:local_solve] == :Optimal || m.status[:local_solve] == :Suboptimal || m.status[:local_solve] == :UserLimit
        # bound tightening goes here - done if local solve is feasible - requires only obj value
    else
        error("NLP local solve is $(m.status[:local_solve]) - quitting solve.")
        quit()
        # do bound tightening without objective value
        # local_solve(m) to generate a feasible solution which is a starting point for bounding_solve
        # if this does not produce an feasible solution then solve atmc without discretization and use as a starting point
    end

    initialize_discretization(m)
    cputime_presolve = time() - start_presolve
    m.logs[:presolve_time] += cputime_presolve
    (m.log_level > 0) && println("Presolve ended.")
    (m.log_level > 0) && println("Presolve time = $(round(m.logs[:total_time],2))s")
    return

end

#=
    Main Adaptive Partitioning (AP) Algorithm.
=#
function global_solve(m::PODNonlinearModel; kwargs...)

    logging_head()
    while m.best_rel_gap > m.rel_gap && m.logs[:time_left] > 0.0
        m.logs[:n_iter] += 1
        update_time_limit(m)        # Updates the time limit, if that option is supplied, TILIM = TILIM - m.logs[:total_time]
        create_bounding_mip(m)      # Build the bounding ATMC model
        bounding_solve(m)           # Solve bounding model
        add_discretization(m)       # Add extra discretizations
        local_solve(m)              # Solve upper bounding model
        m.best_rel_gap = (m.best_obj - m.best_bound)/m.best_obj
        if m.log_level > 0
            logging_row_entry(m)
        end
        if m.logs[:n_iter] > 10
            error("TESTING...")
        end
    end

end

function local_solve(m::PODNonlinearModel; presolve = false)

    convertor = Dict(:Max=>:>, :Min=>:<)
    local_solve_nlp_model = MathProgBase.NonlinearModel(m.nlp_local_solver)
    if presolve == false
        l_var, u_var = tighten_bounds(m)
    end
    MathProgBase.loadproblem!(local_solve_nlp_model, m.num_var_orig, m.num_constr_orig, m.l_var_orig, m.u_var_orig, m.l_constr_orig, m.u_constr_orig, m.sense_orig, m.d_orig)
    MathProgBase.setwarmstart!(local_solve_nlp_model, m.best_sol[1:m.num_var_orig])

    start_local_solve = time()
    MathProgBase.optimize!(local_solve_nlp_model)
    cputime_local_solve = time() - start_local_solve
    m.logs[:total_time] += cputime_local_solve
    m.logs[:time_left] = max(0.0, m.logs[:time_left] - cputime_local_solve)

    local_solve_nlp_status = MathProgBase.status(local_solve_nlp_model)
    if local_solve_nlp_status == :Optimal || local_solve_nlp_status == :Suboptimal || local_solve_nlp_status == :UserLimit
        candidate_obj = MathProgBase.getobjval(local_solve_nlp_model)
        if eval(convertor[m.sense_orig])(candidate_obj, m.best_obj + 1e-10)
            m.best_obj = candidate_obj
            m.best_sol = MathProgBase.getsolution(local_solve_nlp_model)
            m.sol_incumb_ub = copy(m.best_sol)  # temp holder can be removed
            m.status[:feasible_solution] = :Detected
        end
        m.status[:local_solve] = local_solve_nlp_status
        return
    elseif local_solve_nlp_status == :Infeasible
        (m.log_level > 0) && (presolve == true) && println("[PRESOLVE] NLP local solve is infeasible.")
        (m.log_level > 0) && (presolve == false) && println("Incumbent unchanged due to infeasible local solve.")
        m.status[:local_solve] = :Infeasible
        return
    elseif local_solve_nlp_status == :Unbounded
        (presolve == true) && warn("[PRESOLVE] NLP local solve is unbounded.")
        (presolve == false) && warn("[LOCAL SOLVE] NLP local solve is unbounded.")
        m.status[:local_solve] = :Unbounded
        return
    else
        (presolve == true) && error("[PRESOLVE] NLP local solve failure.")
        (presolve == false) && warn("[LOCAL SOLVE] NLP local solve failure.")
        m.status[:local_solve] = :Error
        return
    end

    return
end

function bounding_solve(m::PODNonlinearModel; kwargs...)

    convertor = Dict(:Max=>:<, :Min=>:>)
    update_time_limit(m)
    start_bounding_solve = time()
    status = solve(m.model_mip, suppress_warnings=true)
    cputime_bounding_solve = time() - start_bounding_solve
    m.logs[:total_time] += cputime_bounding_solve

    if status == :Optimal || status == :Suboptimal || status == :UserLimit
        candidate_bound = getobjectivevalue(m.model_mip)
        if eval(convertor[m.sense_orig])(candidate_bound, m.best_bound + 1e-10)
            m.best_bound = candidate_bound
            m.best_bound_sol = [getvalue(Variable(m.model_mip, i)) for i in 1:m.num_var_orig]
            m.sol_incumb_lb = [getvalue(Variable(m.model_mip, i)) for i in 1:m.num_var_orig] # can remove this
            m.status[:bounding_solve] = status
            m.status[:bound] = :Detected
        end
    elseif status == :Infeasible
        # print_iis_gurobi(m.model_mip)  # Used for debugging
        #= Case when this happens::
            1) :UserLimits ? but will this return :UserLimits?
            2) :Numerical Difficult
        =#
        m.status[:bounding_solve] = status
        error("[MIP INFEASIBLE] There is some issue about LB problem")
    elseif status == :Unbounded
        m.status[:bounding_solve] = status
        error("[MIP UNBOUNDED] MIP solver failure")
    else
        error("[MIP UNEXPECTED] MIP solver failure.")
    end
end

#=========================================================
 Logging and printing functions
=========================================================#

# Create dictionary of logs for timing and iteration counts
function create_logs!(m)

    logs = Dict{Symbol,Any}()

    # Timers
    logs[:presolve_time] = 0.       # Total presolve-time of the algorithm
    logs[:total_time] = 0.          # Total run-time of the algorithm
    logs[:time_left] = m.timeout    # Total remaining time of the algorithm if timeout is specified

    # Counters
    logs[:n_iter] = 0           # Number of iterations in iterative
    logs[:n_feas] = 0           # Number of times get a new feasible solution
    logs[:ub_incumb_cnt] = 0    # Number of incumbents detected on upper bound
    logs[:lb_incumb_cnt] = 0    # Number of incumebnts detected on lower bound

    m.logs = logs
end

function logging_summary(m::PODNonlinearModel)
    if m.log_level > 0
        @printf "full problem loaded into POD.\n"
        @printf "number of constraints = %d.\n" m.num_constr_orig
        @printf "number of non-linear constraints = %d.\n" m.num_nlconstr_orig
        @printf "number of linear constraints = %d.\n" m.num_lconstr_orig
        @printf "number of variables = %d.\n" m.num_var_orig

        @printf "relative optimality gap criteria = %.5f (%.4f %%)\n" m.rel_gap (m.rel_gap*100)
        @printf "algorithm for selecting variables to discretize = %d\n" m.var_discretization_algo
    end
end

function logging_head()
    println(" | FS        | BD        | GAP\%      | USED TIME | TIME LEFT | Iter ")
end

function logging_row_entry(m::PODNonlinearModel; kwargs...)
    b_len = 10
    UB_block = string(" ", round(m.best_obj,4), " " ^ (b_len - length(string(round(m.best_obj, 4)))))
    LB_block = string(" ", round(m.best_bound,4), " " ^ (b_len - length(string(round(m.best_bound, 4)))))
    GAP_block = string(" ", round(m.best_rel_gap*100,5), " " ^ (b_len - length(string(round(m.best_rel_gap*100,5)))))
    UTIME_block = string(" ", round(m.logs[:total_time],2), "s", " " ^ (b_len - 1 - length(string(round(m.logs[:total_time],2)))))
    LTIME_block = string(" ", round(m.logs[:time_left],2), "s", " " ^ (b_len - 1 - length(string(round(m.logs[:time_left],2)))))
    ITER_block = string(" ", m.logs[:n_iter])
    println(" |",UB_block,"|",LB_block,"|",GAP_block,"|",UTIME_block,"|",LTIME_block,"|",ITER_block)
end

# Create dictionary of statuses for POD algorithm
function create_status!(m)

    status = Dict{Symbol,Symbol}()

    status[:presolve] = :none                   # Status of presolve
    status[:local_solve] = :none                # Status of local solve
    status[:bounding_solve] = :none              # Status of bounding solve
    status[:lower_bounding_solve] = :none        # Status of lower bonding solve
    status[:upper_bounding_solve] = :none       # Status of bounding solve
    status[:feasible_solution] = :none          # Status of whether a upper bound is detected or not
    status[:upper_bound] = :none                # Status of whether a upper bound has been detected
    status[:lower_bound] = :none                # Status of whether a lower bound has been detected
    status[:bound] = :none                      # Status of whether a bound has been detected
    status[:bound_tightening_solve] = :none    # Status of bound-tightening solve

    m.status = status
end

function summary_status(m::PODNonlinearModel)

    if m.status[:bound] == :Detected && m.status[:feasible_solution] == :Detected
        if m.best_rel_gap >= m.rel_gap
            m.pod_status = :UserLimits
        else
            m.pod_status = :Optimal
        end
    elseif m.status[:bound] == :Detected && m.status[:feasible_solution] == :none
        m.pod_status = :Infeasible
    elseif m.status[:bound] == :none && m.status[:feasible_solution] == :Detected
        m.pod_status = :Heuristic
    else
        error("[UNEXPECTED] Missing bound and feasible solution during status summary.")
    end

end
