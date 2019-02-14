

function create_default_solver_options()
    solver_options = Dict{Symbol, Union{Nothing, MOI.AbstractOptimizer, Any}}(
        :nlp_optimizer => nothing,          
        :minlp_optimizer => nothing,          
        :mip_optimizer => nothing,       
        :log_level => 1,                                              
        :time_limit => Inf, 
        :max_iter => 100,
        :rel_gap => 1e-4, 
        :abs_gap => 1e-6, 
        :opt_tol => 1e-6,
        :disc_partition_method => :adaptive,
        :disc_ratio => 4.0,
        :disc_uniform_rate => 2,
        :disc_divert_chunks => 5,
        :disc_abs_width_tol => 1e-4,
        :disc_rel_width_tol => 1e-4,
        :disc_consecutive_forbid => 0, 
        :bt => true, 
        :presolve_time_limit => 900.0,
        :bt_max_iter => 10,
        :bt_width_tol => 1e-3,
        :bt_output_tol => 1e-5,
        :bt_algo => :continuous,
        :bt_relax => false, 
        :bt_mip_time_limit => 200.0, 
        :bp => true,
        :is_convex => false 
    )
    return solver_options 
end 

function update_solver_options(default_options, given_options) 
    available_options = keys(default_options)
    for (key, val) in given_options 
        if !(key in available_options)
            error(LOGGER, "the option $key is not supported by Alpine.Optimizer; check the documentation for supported options")
        else 
            default_options[key] = val 
        end 
    end 
    return 
end 
