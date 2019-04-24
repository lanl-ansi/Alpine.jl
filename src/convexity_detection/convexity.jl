
""" 
    run_convexity_detection!(model::MOI.AbstractOptimizer)

Entry point for the convexity detection routines  
"""

function run_convexity_detection!(model::MOI.AbstractOptimizer)
    run_quadratic_convexity_detection!(model)
    run_quadratic_nl_convexity_detection!(model)
    # dag_convexity_propagation!(model)
end 

