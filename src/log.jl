# Create dictionary of logs for timing and iteration counts
function create_logs!(m)

    logs = Dict{Symbol,Any}()

    # Timers
    logs[:presolve_time] = 0.       # Total presolve-time of the algorithm
    logs[:total_time] = 0.          # Total run-time of the algorithm
    logs[:time_left] = m.timeout    # Total remaining time of the algorithm if timeout is specified

    # Values
    logs[:obj] = []                 # Iteration based objective
    logs[:bound] = []               # Iteration based objective

    # Counters
    logs[:n_iter] = 0               # Number of iterations in iterative
    logs[:n_feas] = 0               # Number of times get a new feasible solution
    logs[:ub_incumb_cnt] = 0        # Number of incumbents detected on upper bound
    logs[:lb_incumb_cnt] = 0        # Number of incumebnts detected on lower bound
    logs[:bt_iter] = 0

    m.logs = logs
end

function reset_timer(m::AlpineNonlinearModel)
    m.logs[:total_time] = 0.
    m.logs[:time_left] = m.timeout
    return m
end

function logging_summary(m::AlpineNonlinearModel)

    if m.loglevel > 0
        printstyled(:light_yellow, "full problem loaded into Alpine\n")
        println("problen sense $(m.sense_orig)")
        println("# of constraints = ", m.num_constr_orig)
        println("# of non-linear constraints = ", m.num_nlconstr_orig)
        println("# of linear constraints = ", m.num_lconstr_orig)
        println("# of continuous variables = ", length([i for i in 1:m.num_var_orig if m.var_type[i] == :Cont]))
        println("# of binary variables = ", length([i for i in 1:m.num_var_orig if m.var_type[i] == :Bin]))
        println("# of integer variables = ", length([i for i in 1:m.num_var_orig if m.var_type[i] == :Int]))

        println("detected nonlinear terms = ", length(m.nonconvex_terms))
        for i in ALPINE_C_NLTERMS
            cnt = length([1 for j in keys(m.nonconvex_terms) if m.nonconvex_terms[j][:nonlinear_type] == i])
            cnt > 0 && println("\tTerm $(i) Count = $(cnt) ")
        end
        println("# of variables involved in nonlinear terms = ", length(m.candidate_disc_vars))
        println("# of selected variables to discretize = ", length(m.disc_vars))

        println("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *")
        println("                        SUB-SOLVERS  ")
        m.minlp_solver != UnsetSolver() && println("MINLP local solver = ", split(string(m.minlp_solver),".")[1])
        println("NLP local solver = ", split(string(m.nlp_solver),".")[1])
        println("MIP solver = ", split(string(m.mip_solver),".")[1])
        println("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *")
        println("                    SOLVER CONFIGURATION ")
        println("maximum solution time = ", m.timeout)
        println("maximum iterations =  ", m.maxiter)
        @printf "relative optimality gap criteria = %.5f (%.4f %%)\n" m.relgap (m.relgap*100)
        println("default tolerance = ", m.tol)
        m.recognize_convex && println("actively recognize convex patterns")
        m.recognize_convex && println("# of detected convex constraints = $(length([i for i in m.constr_structure if i == :convex]))")
        println("basic bound propagation = ", m.presolve_bp)
        println("use piece-wise relaxation formulation on integer variables = ", m.int_enable)
        println("piece-wise relaxation formulation = $(m.convhull_formulation) formulation")
        println("method for picking discretization variable = $(m.disc_var_pick)")
        println("conseuctive solution rejection = after ", m.disc_consecutive_forbid, " times")
        if m.disc_ratio_branch
            println("discretization ratio branch activated")
        else
            println("discretization ratio = ", m.disc_ratio)
        end
        (m.convhull_ebd) && println("using convhull_ebd formulation")
        (m.convhull_ebd) && println("encoding method = $(m.convhull_ebd_encode)")
        (m.convhull_ebd) && println("independent branching scheme = $(m.convhull_ebd_ibs)")
        println("bound tightening presolver = ", m.presolve_bt)
        m.presolve_bt && println("bound tightening presolve maximum iteration = ", m.presolve_maxiter)
        m.presolve_bt && println("bound tightening presolve algorithm = ", m.presolve_bt_algo)
        m.presolve_bt && println("bound tightening presolve width tolerance = ", m.presolve_bt_width_tol)
        m.presolve_bt && println("bound tightening presolve output tolerance = ", m.presolve_bt_output_tol)
        m.presolve_bt && println("bound tightening presolve relaxation = ", m.presolve_bt_relax)
        m.presolve_bt && println("bound tightening presolve mip regulation time = ", m.presolve_bt_mip_timeout)
        println("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *")
    end

    # Additional warnings
    m.mip_solver_id == "Gurobi" && @warn "Alpine support Gurobi solver 7.0+ ..."
end

function logging_head(m::AlpineNonlinearModel)
	if m.logs[:time_left] < Inf
		printstyled(:light_yellow, " | NLP           | MIP           || Objective     | Bound         | GAP (%)       | CLOCK         | TIME LEFT     | Iter   \n")
	else
		printstyled(:light_yellow, " | NLP           | MIP           || Objective     | Bound         | GAP (%)       | CLOCK         | | Iter   \n")
	end
end

function logging_row_entry(m::AlpineNonlinearModel; kwargs...)

    options = Dict(kwargs)

    b_len = 14
    if !isempty(m.logs[:obj]) && isa(m.logs[:obj][end], Float64)
        objstr = string(round(m.logs[:obj][end]; digits=4))
        spc = max(0, b_len - length(objstr))
    else
        objstr = string("-")
        spc = max(0, b_len - length(objstr))
    end
    UB_block = string(" ", objstr, " " ^ spc)

    if isa(m.logs[:bound][end], Float64)
        bdstr = string(round(m.logs[:bound][end]; digits=4))
        spc = max(0, b_len - length(bdstr))
    else
        bdstr = string(m.logs[:bound][end])
        spc = b_len - length(bdstr)
    end
    LB_block = string(" ", bdstr, " " ^ spc)

    bobjstr = string(round(m.best_obj; digits=4))
    spc = max(0, b_len - length(bobjstr))
    incumb_UB_block = string(" ", bobjstr, " " ^ spc)

    bbdstr = string(round(m.best_bound; digits=4))
    spc = max(0, b_len - length(bbdstr))
    incumb_LB_block = string(" ", bbdstr , " " ^ spc)

    rel_gap = round(m.best_rel_gap*100, digits=5)
    rel_gap > 999 ? rel_gap = "LARGE" : rel_gap = string(rel_gap)
    GAP_block = string(" ", rel_gap, " " ^ (b_len - length(rel_gap)))

    UTIME_block = string(" ", round(m.logs[:total_time]; digits=2), "s", " " ^ (b_len - 1 - length(string(round(m.logs[:total_time]; digits=2)))))

    if m.logs[:time_left] < Inf
		LTIME_block = string(" ", round(m.logs[:time_left]; digits=2), "s", " " ^ (b_len - 1 - length(string(round(m.logs[:time_left]; digits=2)))))
	else
		LTIME_block = " "
	end

    haskey(options, :finsih_entry) ? (ITER_block = string(" ", "finish")) : (ITER_block = string(" ", m.logs[:n_iter]))

    if m.colorful_pod == "random"
        colors = [:blue, :cyan, :green, :red, :light_red, :light_blue, :light_cyan, :light_green, :light_magenta, :light_re, :light_yellow, :white, :yellow]
        print(" |")
        printstyled(rand(colors),UB_block)
        print("|")
        printstyled(rand(colors),LB_block)
        print("||")
        printstyled(rand(colors),incumb_UB_block)
        print("|")
        printstyled(rand(colors),incumb_LB_block)
        print("|")
        printstyled(rand(colors),GAP_block)
        print("|")
        printstyled(rand(colors),UTIME_block)
        print("|")
        printstyled(rand(colors),LTIME_block)
        print("|")
        printstyled(rand(colors),ITER_block)
        print("\n")
    elseif m.colorful_pod == "solarized"
        printstyled(m.logs[:n_iter]+21, " |$(UB_block)|$(LB_block)||$(incumb_UB_block)|$(incumb_LB_block)|$(GAP_block)|$(UTIME_block)|$(LTIME_block)|$(ITER_block)\n")
    elseif m.colorful_pod == "warmer"
        printstyled(max(20,170-m.logs[:n_iter]), " |$(UB_block)|$(LB_block)||$(incumb_UB_block)|$(incumb_LB_block)|$(GAP_block)|$(UTIME_block)|$(LTIME_block)|$(ITER_block)\n")
    elseif m.colorful_pod == false
        println(" |",UB_block,"|",LB_block,"||",incumb_UB_block,"|",incumb_LB_block,"|",GAP_block,"|",UTIME_block,"|",LTIME_block,"|",ITER_block)
    end

    return
end


#=========================================================
 Logging and printing functions
=========================================================#

# Create dictionary of statuses for Alpine algorithm
function create_status!(m)

    status = Dict{Symbol,Symbol}()

    status[:presolve] = :none                   # Status of presolve
    status[:local_solve] = :none                # Status of local solve
    status[:bounding_solve] = :none             # Status of bounding solve
    status[:feasible_solution] = :none          # Status of whether a upper bound is detected or not
    status[:bound] = :none                      # Status of whether a bound has been detected

    m.status = status
end


"""
    This function summarize the eventual solver status based on all available infomration
        recorded in the solver. The output status is self-defined which requires users to
        read our documentation to understand what details is behind a status symbol.
"""
function summary_status(m::AlpineNonlinearModel)

    # Alpine Solver Status Definition
    # :Optimal : normal termination with gap closed within time limits
    # :UserLimits : any non-optimal termination related to user-defined parameters
    # :Infeasible : termination with relaxation proven infeasible or detection of
    #               variable bound conflicts
    # :Heuristic : termination with feasible solution found but not bounds detected
    #               happens when lower bound problem is extremely hard to solve
    # :Unknown : termination with no exception recorded

    if m.status[:bound] == :Detected && m.status[:feasible_solution] == :Detected
        m.best_rel_gap > m.relgap ? m.pod_status = :UserLimits : m.pod_status = :Optimal
    elseif m.status[:bounding_solve] == :Infeasible
        m.pod_status = :Infeasible
    elseif m.status[:bound] == :Detected && m.status[:feasible_solution] == :none
        m.pod_status = :UserLimits
    elseif m.status[:bound] == :none && m.status[:feasible_solution] == :Detected
        m.pod_status = :Heuristic
    else
        @warn "[EXCEPTION] Indefinite Alpine status. Please report your instance and configuration as and Issue (https://github.com/lanl-ansi/Alpine.jl/issues) to help us make Alpine better."
    end

    printstyled(:green, "\n Alpine ended with status $(m.pod_status)\n")

    return
end
