# Create dictionary of logs for timing and iteration counts
function create_logs!(m)
    logs = Dict{Symbol,Any}()

    # Timers
    logs[:presolve_time] = 0.0                           # Total presolve-time of the algorithm
    logs[:total_time] = 0.0                              # Total run-time of the algorithm
    logs[:time_left] = get_option(m, :time_limit)    # Total remaining time of the algorithm if time-out is specified

    # Values
    logs[:obj] = []                 # Iteration-based objective
    logs[:bound] = []               # Iteration-based objective bound

    # Counters
    logs[:n_iter] = 0               # Number of iterations
    logs[:n_feas] = 0               # Number of times a new feasible solution is obtained
    logs[:ub_incumb_cnt] = 0        # Number of incumbents detected in the upper bound
    logs[:lb_incumb_cnt] = 0        # Number of incumebnts detected in the lower bound
    logs[:bt_iter] = 0

    return m.logs = logs
end

function reset_timer(m::Optimizer)
    m.logs[:total_time] = 0.0
    m.logs[:time_left] = get_option(m, :time_limit)
    return m
end

function logging_summary(m::Optimizer)
    if get_option(m, :log_level) > 0
        printstyled("\nPROBLEM STATISTICS\n", color = :cyan, bold = true)
        is_min_sense(m) && (println("  Objective sense = Min"))
        is_max_sense(m) && (println("  Objective sense = Max"))
        println(
            "  # Variables = ",
            length([i for i in 1:m.num_var_orig if m.var_type[i] == :Cont]) +
            length([i for i in 1:m.num_var_orig if m.var_type[i] == :Bin]) +
            length([i for i in 1:m.num_var_orig if m.var_type[i] == :Int]),
        )
        println(
            "  # Bin-Int Variables = ",
            length([i for i in 1:m.num_var_orig if m.var_type[i] == :Bin]) + length([i for i in 1:m.num_var_orig if m.var_type[i] == :Int]),
        )
        println("  # Constraints = ", m.num_constr_orig)
        println("  # NL Constraints = ", m.num_nlconstr_orig)
        println("  # Linear Constraints = ", m.num_lconstr_orig)
        get_option(m, :recognize_convex) && println(
            "  # Detected convex constraints = $(length([i for i in m.constr_structure if i == :convex]))",
        )
        println("  # Detected nonlinear terms = ", length(m.nonconvex_terms))
        println(
            "  # Variables involved in nonlinear terms = ",
            length(m.candidate_disc_vars),
        )
        println("  # Potential variables for partitioning = ", length(m.disc_vars))

        printstyled("SUB-SOLVERS USED BY ALPINE\n", color = :cyan, bold = true)
        if get_option(m, :minlp_solver) === nothing
            println("  NLP local solver = ", m.nlp_solver_id)
        else
            println("  MINLP local solver = ", m.minlp_solver_id)
        end
        println("  MIP solver = ", m.mip_solver_id)

        printstyled("ALPINE CONFIGURATION\n", color = :cyan, bold = true)
        project = read(joinpath(dirname(@__DIR__), "Project.toml"), String)
        m_version = match(r"version \= \"(.+?)\"", project)
        println("  Alpine version = ", m_version[1])
        if is_min_sense(m)
            println(
                "  Maximum iterations (lower-bounding MIPs) = ",
                get_option(m, :max_iter),
            )
        elseif is_max_sense(m)
            println(
                "  Maximum iterations (upper-bounding MIPs) = ",
                get_option(m, :max_iter),
            )
        else
            println("  Maximum iterations (bounding MIPs) =  ", get_option(m, :max_iter))
        end

        println("  Relative global optimality gap = ", get_option(m, :rel_gap) * 100, "%")

        if get_option(m, :disc_var_pick) == 0
            println("  Potential variables chosen for partitioning = All")
        elseif get_option(m, :disc_var_pick) == 1
            println(
                "  Potential variables chosen for partitioning = Minimum vertex cover",
            )
        end

        if get_option(m, :partition_scaling_factor_branch)
            println("  Partition scaling factor branch activated")
        else
            println(
                "  Partition scaling factor = ",
                get_option(m, :partition_scaling_factor),
            )
        end
        (get_option(m, :convhull_ebd)) && println("  Using convhull_ebd formulation")
        (get_option(m, :convhull_ebd)) &&
            println("  Encoding method = $(get_option(m, :convhull_ebd_encode))")
        (get_option(m, :convhull_ebd)) && println(
            "  Independent branching scheme = $(get_option(m, :convhull_ebd_ibs))",
        )
        println("  Bound-tightening presolve = ", get_option(m, :presolve_bt))
        get_option(m, :presolve_bt) && println(
            "  Maximum iterations (OBBT) = ",
            get_option(m, :presolve_bt_max_iter),
        )
    end
end

function logging_head(m::Optimizer)
    if is_min_sense(m)
        printstyled("LOWER-BOUNDING ITERATIONS", color = :cyan, bold = true)
        UB_iter = "Incumbent"
        UB = "Best Incumbent"
        LB = "Lower Bound"
    elseif is_max_sense(m)
        printstyled("UPPER-BOUNDING ITERATIONS", color = :cyan, bold = true)
        UB_iter = "Incumbent"
        UB = "Best Incumbent"
        LB = "Upper Bound"
    end
    println(
        "\n====================================================================================================",
    )
    if m.logs[:time_left] < Inf
        printstyled(
            bold = true,
            "| Iter   | $UB_iter       | $UB      | $LB        | Gap (%)         | Time      \n",
        )
    end
end

function logging_row_entry(m::Optimizer; kwargs...)
    options = Dict(kwargs)

    b_len = 16
    if !isempty(m.logs[:obj]) && isa(m.logs[:obj][end], Float64)
        objstr = string(round(m.logs[:obj][end]; digits = 4))
        spc = max(0, b_len - length(objstr))
    else
        objstr = string("-")
        spc = max(0, b_len - length(objstr))
    end
    UB_block = string(" ", objstr, " "^spc)

    if expr_isconst(m.obj_expr_orig)
        bdstr = eval(m.obj_expr_orig)
        spc = b_len - length(bdstr)
    elseif isa(m.logs[:bound][end], Float64)
        bdstr = string(round(m.logs[:bound][end]; digits = 4))
        spc = max(0, b_len - length(bdstr))
    else
        bdstr = string(m.logs[:bound][end])
        spc = b_len - length(bdstr)
    end
    LB_block = string(" ", bdstr, " "^spc)

    bobjstr = string(round(m.best_obj; digits = 4))
    spc = max(0, b_len + 4 - length(bobjstr))
    incumb_UB_block = string(" ", bobjstr, " "^spc)

    bbdstr = string(round(m.best_bound; digits = 4))
    spc = max(0, b_len + 3 - length(bbdstr))
    incumb_LB_block = string(" ", bbdstr, " "^spc)

    rel_gap = round(m.best_rel_gap * 100, digits = 3)
    rel_gap > 999 ? rel_gap = "LARGE" : rel_gap = string(rel_gap)
    GAP_block = string(" ", rel_gap, " "^(b_len - length(rel_gap)))

    UTIME_block = string(
        " ",
        round(m.logs[:total_time]; digits = 2),
        "s",
        " "^(b_len - 1 - length(string(round(m.logs[:total_time]; digits = 2)))),
    )

    if m.logs[:time_left] < Inf
        LTIME_block = " "
    end

    haskey(options, :finish_entry) ? (ITER_block = string(" ", "finish ")) :
    (ITER_block = string(" ", m.logs[:n_iter], " "^(7 - length(string(m.logs[:n_iter])))))

    println(
        "|",
        ITER_block,
        "|",
        UB_block,
        "|",
        incumb_UB_block,
        "|",
        incumb_LB_block,
        "|",
        GAP_block,
        "|",
        UTIME_block,
        LTIME_block,
    )
    return
end

#Logging and printing functions

# Create dictionary of statuses for Alpine algorithm
function create_status!(m)
    status = Dict{Symbol,MOI.TerminationStatusCode}()

    status[:local_solve] = MOI.OPTIMIZE_NOT_CALLED # Status of local solve
    status[:bounding_solve] = MOI.OPTIMIZE_NOT_CALLED # Status of bounding solve
    m.detected_incumbent = false
    m.detected_bound = false

    return m.status = status
end

"""
This function summarizes the eventual solver status based on all available information
recorded in the solver. The output status is self-defined which requires users to
read our documentation to understand the details behind every status symbols.
"""
function summary_status(m::Optimizer)

    # Alpine Solver Status Definition
    # :Optimal : normal termination with optimality gap closed within time limits
    # :UserLimits : any non-optimal termination related to user-defined parameters
    # :Infeasible : termination with relaxation proven infeasible or detection of
    #               variable bound conflicts
    # :Heuristic : termination with feasible solution found but not bounds detected
    #               happens when lower bound problem is extremely hard to solve
    # :Unknown : termination with no exception recorded

    if m.detected_bound && m.detected_incumbent
        m.alpine_status =
            m.best_rel_gap > get_option(m, :rel_gap) ? MOI.OTHER_LIMIT : MOI.OPTIMAL
    elseif m.status[:bounding_solve] == MOI.INFEASIBLE
        m.alpine_status = MOI.INFEASIBLE
    elseif m.detected_bound && !m.detected_incumbent
        m.alpine_status = MOI.OTHER_LIMIT
    elseif !m.detected_bound && m.detected_incumbent
        m.alpine_status = MOI.LOCALLY_SOLVED
    else
        @warn "  [EXCEPTION] Indefinite Alpine status. Please report your instance (& solver configuration) as an issue (https://github.com/lanl-ansi/Alpine.jl/issues) to help us make Alpine better."
    end

    raw = MOI.get(m, MOI.RawStatusString())
    printstyled("\n*** $raw ***\n")

    return
end

# Some useful logging details:

# for i in ALPINE_C_NLTERMS
#     cnt = length([1 for j in keys(m.nonconvex_terms) if m.nonconvex_terms[j][:nonlinear_type] == i])
#     cnt > 0 && println("\tTerm $(i) Count = $(cnt) ")
# end
# println("  Maximum solution time = ", get_option(m, :time_limit))
# println("  Basic bound propagation = ", get_option(m, :presolve_bp))
# println("  Conseuctive solution rejection = after ", get_option(m, :disc_consecutive_forbid), " times")
# get_option(m, :presolve_bt) && println("bound tightening presolve algorithm = ", get_option(m, :presolve_bt)_algo)
# get_option(m, :presolve_bt) && println("bound tightening presolve width tolerance = ", get_option(m, :presolve_bt)_width_tol)
# get_option(m, :presolve_bt) && println("bound tightening presolve output tolerance = ", get_option(m, :presolve_bt)_output_tol)
# get_option(m, :presolve_bt) && println("bound tightening presolve relaxation = ", get_option(m, :presolve_bt)_relax)
# get_option(m, :presolve_bt) && println("bound tightening presolve mip regulation time = ", get_option(m, :presolve_bt)_mip_time_limit)
