

function get_default_solver_options()
    nlp_optimizer = nothing 
    minlp_optimizer = nothing 
    mip_optimizer = nothing 
    log_level = 1 
    time_limit = Inf 
    max_iter = 100 
    rel_gap = 1e-4 
    abs_gap = 1e-6 
    opt_tol = 1e-6 
    disc_partition_method = :adaptive 
    disc_ratio = 4.0 
    disc_uniform_rate = 2
    disc_divert_chunks = 5 
    disc_abs_width_tol = 1e-4 
    disc_rel_width_tol = 1e-4
    disc_consecutive_forbid = 0
    bt = true 
    presolve_time_limit = 900.0
    bt_max_iter = 10
    bt_width_tol = 1e-3
    bt_improvement_tol = 1e-2
    bt_precision = 5
    bt_algo = :continuous
    bt_relax = false 
    bt_mip_time_limit = 200.0 
    bp = true
    bp_max_iter = 3
    is_problem_convex = false 
    perform_bt_only = false
    perform_bounding_solve_only = false
    
    return SolverOptions(
        nlp_optimizer, minlp_optimizer, mip_optimizer, 
        log_level, time_limit, max_iter, 
        rel_gap, abs_gap, opt_tol, 
        disc_partition_method, disc_ratio, 
        disc_uniform_rate, disc_divert_chunks, 
        disc_abs_width_tol, disc_rel_width_tol, 
        disc_consecutive_forbid, 
        bt, presolve_time_limit, 
        bt_max_iter, bt_width_tol, bt_improvement_tol,
        bt_precision, bt_algo, bt_relax, bt_mip_time_limit, 
        bp, bp_max_iter, is_problem_convex, 
        perform_bt_only, perform_bounding_solve_only
    )
end 


function combine_options(options)
    options_dict = Dict{Symbol,Any}()
    for kv in options
        if !in(kv[1], fieldnames(SolverOptions))
            warn(LOGGER, "Option " * string(kv[1]) * " is not available")
        end
        options_dict[kv[1]] = kv[2]
    end

    if !haskey(options_dict, :nlp_optimizer) && !haskey(options_dict, :minlp_optimizer)
        error(LOGGER, "either an NLP or a MINLP solver has to be specified using nlp_optimizer or minlp_optimizer option")
    end

    if !haskey(options_dict, :mip_optimizer) 
        error(LOGGER, "MIP solver has to be specified using mip_optimizer option")
    end 

    defaults = get_default_solver_options()
    for fname in fieldnames(SolverOptions)
        if haskey(options_dict, fname)

            if fieldtype(SolverOptions, fname) != typeof(options_dict[fname])
                options_dict[fname] = convert(fieldtype(SolverOptions, fname), options_dict[fname])
            end
            setfield!(defaults, fname, options_dict[fname])
        end
    end
    return defaults
end 

