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

function reset_timer(m::PODNonlinearModel)
    m.logs[:total_time] = 0.
    m.logs[:time_left] = m.timeout
    return m
end

function logging_summary(m::PODNonlinearModel)

    if m.loglevel > 0
        print_with_color(:light_yellow, "full problem loaded into POD\n")
        println("problen sense $(m.sense_orig)")
        @printf "number of constraints = %d\n" m.num_constr_orig
        @printf "number of non-linear constraints = %d\n" m.num_nlconstr_orig
        @printf "number of linear constraints = %d\n" m.num_lconstr_orig
        @printf "number of variables = %d\n" m.num_var_orig

        m.minlp_solver != UnsetSolver() && println("MINLP local solver = ", split(string(m.minlp_solver),".")[1])
        println("NLP local solver = ", split(string(m.nlp_solver),".")[1])
        println("MIP solver = ", split(string(m.mip_solver),".")[1])

        println("maximum solution time = ", m.timeout)
        println("maximum iterations =  ", m.maxiter)
        @printf "relative optimality gap criteria = %.5f (%.4f %%)\n" m.relgap (m.relgap*100)
        println("detected nonlinear terms = $(length(m.nonlinear_terms))")
        println("number of variables involved in nonlinear terms = $(length(m.all_nonlinear_vars))")
        println("number of selected variables to discretize = $(length(m.var_disc_mip))")

        for i in POD_C_NLTERMS
            cnt = length([1 for j in keys(m.nonlinear_terms) if m.nonlinear_terms[j][:nonlinear_type] == i])
            cnt > 0 && println("Term $(i) Count = $(cnt) ")
        end

        m.bilinear_convexhull && println("bilinear treatment = convex hull formulation")
        m.monomial_convexhull && println("monomial treatment = convex hull formulation")

        println("using piece-wise convex hull : $(m.convhull_formulation) formulation")
        println("using method $(m.disc_var_pick) for picking discretization variable...")

        (m.convhull_ebd) && println("using convhull_ebd formulation")
        (m.convhull_ebd) && println("encoding method = $(m.convhull_ebd_encode)")
        (m.convhull_ebd) && println("independent branching scheme = $(m.convhull_ebd_ibs)")
    end

    # Additional warnings
    m.mip_solver_id == "Gurobi" && warn("POD support Gurobi solver 7.0+ ...")
end

function logging_head(m::PODNonlinearModel)
	if m.logs[:time_left] < Inf
		print_with_color(:light_yellow, " | NLP           | MIP           || Objective     | Bound         | GAP\%          | CLOCK         | TIME LEFT     | Iter   \n")
	else
		print_with_color(:light_yellow, " | NLP           | MIP           || Objective     | Bound         | GAP\%          | CLOCK         | | Iter   \n")
	end

	return
end

function logging_row_entry(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    b_len = 14
    if isa(m.logs[:obj][end], Float64)
        UB_block = string(" ", round(m.logs[:obj][end],4), " " ^ (b_len - length(string(round(m.logs[:obj][end], 4)))))
    else
        UB_block = string(" ", string(m.logs[:obj][end]), " " ^ (b_len - length(string(m.logs[:obj][end]))))
    end
    LB_block = string(" ", round(m.logs[:bound][end],4), " " ^ (b_len - length(string(round(m.logs[:bound][end], 4)))))
    incumb_UB_block = string(" ", round(m.best_obj,4), " " ^ (b_len - length(string(round(m.best_obj, 4)))))
    incumb_LB_block = string(" ", round(m.best_bound,4), " " ^ (b_len - length(string(round(m.best_bound, 4)))))
    GAP_block = string(" ", round(m.best_rel_gap*100,5), " " ^ (b_len - length(string(round(m.best_rel_gap*100,5)))))
    UTIME_block = string(" ", round(m.logs[:total_time],2), "s", " " ^ (b_len - 1 - length(string(round(m.logs[:total_time],2)))))
    if m.logs[:time_left] < Inf
		LTIME_block = string(" ", round(m.logs[:time_left],2), "s", " " ^ (b_len - 1 - length(string(round(m.logs[:time_left],2)))))
	else
		LTIME_block = " "
	end
    haskey(options, :finsih_entry) ? (ITER_block = string(" ", "finish")) : (ITER_block = string(" ", m.logs[:n_iter]))

    if m.colorful_pod == "random"
        colors = [:blue, :cyan, :green, :red, :light_red, :light_blue, :light_cyan, :light_green, :light_magenta, :light_re, :light_yellow, :white, :yellow]
        print(" |")
        print_with_color(rand(colors),UB_block)
        print("|")
        print_with_color(rand(colors),LB_block)
        print("||")
        print_with_color(rand(colors),incumb_UB_block)
        print("|")
        print_with_color(rand(colors),incumb_LB_block)
        print("|")
        print_with_color(rand(colors),GAP_block)
        print("|")
        print_with_color(rand(colors),UTIME_block)
        print("|")
        print_with_color(rand(colors),LTIME_block)
        print("|")
        print_with_color(rand(colors),ITER_block)
        print("\n")
    elseif m.colorful_pod == "solarized"
        print_with_color(m.logs[:n_iter]+21, " |$(UB_block)|$(LB_block)||$(incumb_UB_block)|$(incumb_LB_block)|$(GAP_block)|$(UTIME_block)|$(LTIME_block)|$(ITER_block)\n")
    elseif m.colorful_pod == "warmer"
        print_with_color(max(20,170-m.logs[:n_iter]), " |$(UB_block)|$(LB_block)||$(incumb_UB_block)|$(incumb_LB_block)|$(GAP_block)|$(UTIME_block)|$(LTIME_block)|$(ITER_block)\n")
    elseif m.colorful_pod == false
        println(" |",UB_block,"|",LB_block,"||",incumb_UB_block,"|",incumb_LB_block,"|",GAP_block,"|",UTIME_block,"|",LTIME_block,"|",ITER_block)
    end

    return
end


#=========================================================
 Logging and printing functions
=========================================================#

# Create dictionary of statuses for POD algorithm
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
function summary_status(m::PODNonlinearModel)

    # POD Solver Status Definition
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
        warn("[EXCEPTION] Indefinite POD status. Please report your instance and configuration as and Issue (https://github.com/lanl-ansi/POD.jl/issues) to help us make POD better.")
    end

    print_with_color(:green, "\n POD ended with status $(m.pod_status)\n")

    return
end
