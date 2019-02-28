"""
Simple convexity detection for quadratic inequality constraints 
"""

function run_convexity_detection!(model::MOI.AbstractOptimizer)
    for constraint in model.quadratic_le_constraints 

    end 

    for constraint in model.quadratic_ge_constraints 

    end

end 
