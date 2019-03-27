"""
Feasibility-based bound-tightening for Alpine.jl
"""
function run_fbbt!(model::MOI.AbstractOptimizer)

    num_iter = model.solver_options.bt_max_iter 
    for i in 1:num_iter 
        run_forward_propagation!(model)
        (is_infeasible(model)) && (return)
        run_backward_propagation!(model)
        (is_infeasible(model)) && (return)
    end 

    return
end 

"""
Forward propagation 
""" 
function run_forward_propagation!(model::MOI.AbstractOptimizer)
    fbbt_forward_linear_constraints!(model)
    # fbbt_forward_quadratic_constraints!(model)
    # fbbt_forward_dag!(model)
    # fbbt_forward_nl_constraints!(model)
    return 
end 

"""
Backward propogation 
"""
function run_backward_propagation!(model::MOI.AbstractOptimizer)
    fbbt_backward_linear_constraints!(model)
    # fbbt_backward_quadratic_constraints!(model)
    # fbbt_backward_dag!(model)
    # fbbt_backward_nl_constraints!(model)
    
    return 
end 