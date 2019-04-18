"""
    presolve!(model::MOI.AbstractOptimizer)

This function is the entry point for the presolver of Alpine.jl. 
The different steps involved in the presolve part of the algorithm
(in order) are as follows: 

1. Perform a local solve of the original model. If the original 
problem is an NLP, then the vanilla algorithm to solve the problem
to local optimality is performed with a multi-start option. If it 
is an MINLP and an ``minlp_optimizer`` is provided as a user option, 
the problem is directly provided to the solver to obtain a locally
optimal solution. If the problem is an MINLP and no MINLP solver is 
provided, no initial local solve is performed. 

2. Once step 1 is complete, feasibility-based bound-tightenening 
is performed to obtain tight variable bounds, objective bounds, 
identify redundant constraints, and detect infeasibilities. 

3. If original problem is an NLP, another multi-start local solve 
is performed using the tightened bounds to compute a better local 
optimal solution, if possible. 

4. Run convexity detection algorithms

5. The relaxation is built for performing optimization-based bound
-tightening algorithm

6. The current local optimal objective computed is then used together
with an optimization-based bound-tightening algorithm to compute
tighter variable bounds.

7. The initial root relaxation model is then built and solved 
to obtain the initial bound to the objective function. 

8. Finally, the statistics of the presolver is logged based on 
what the ``log_level`` option. 
"""

function presolve!(model::MOI.AbstractOptimizer)
    
    # perform local solve 
    local_solve!(model, presolve=true)
    
    # perform FBBT 
    if model.solver_options.bp == true 
        run_fbbt!(model)
        
        # check infeasibility and return 
        if ~isa(model.inner.status, Nothing) && model.inner.status.alpine_status == MOI.INFEASIBLE 
            info(LOGGER, "problem infeasibility detected using bound-propagation")
            return
        end 

        # update the variable bounds and perform another local solve 
        update_variable_bounds!(model)
        local_solve!(model, num_standalone_solves=5)
    end

    # if `perform_bp_only` is set to true, the presolver exits
    if model.solver_options.perform_bp_only 
        return
    end

    # convexity detection is peformed
    run_convexity_detection!(model)


    return 
end 