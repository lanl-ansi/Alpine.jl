"""
    print_line(count::Int) 
"""
function print_line(count::Int)
    for i in 1:count
        print("-")
    end
    print("\n") 

    return
end

"""
    print_problem_summary(ap::AlpineProblem)
"""
function print_problem_summary(ap::AlpineProblem)

    printstyled("\nProblem statistics\n", color=:cyan)
    @printf("%-60s : %5d\n", "Number of variables", ap.num_variables)
    
    if ap.num_binary_variables != 0 
        @printf("%-60s : %5d\n", "Number of binary variables", ap.num_binary_variables)
    end
    
    @printf("%-60s : %5d\n", "Number of constraints", ap.num_constraints)
    
    if ap.num_linear_constraints != 0
        @printf("%-60s : %5d\n", "Number of linear constraints", ap.num_linear_constraints )
    end 

    if ap.num_quadratic_constraints != 0 
        @printf("%-60s : %5d\n", "Number of quadratic constraints", ap.num_quadratic_constraints)
    end 

    if ap.num_soc_constraints != 0
        @printf("%-60s : %5d\n", "Number of SOC constraints", ap.num_soc_constraints)
    end 

    if ap.num_rsoc_constraints != 0 
        @printf("%-60s : %5d\n", "Number of RSOC constraints", ap.num_rsoc_constraints)
    end

    if ap.num_nl_constraints != 0
        @printf("%-60s : %5d\n", "Number of nonlinear constraints", ap.num_nl_constraints)
    end 

    println()

    return
end

"""
    print_presolve_summary(model::MOI.AbstractOptimizer)
"""
function print_presolve_summary(model::MOI.AbstractOptimizer)

    printstyled("\nPresolve statistics\n", color=:cyan)

    binary_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_binary))
    integer_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_integer))
    has_discrete_variables = length(binary_variables) > 0 || length(integer_variables) > 0

    if model.inner.incumbent.status == MOI.LOCALLY_SOLVED 
        @printf("%-60s : %10.4f\n", "Initial feasible objective value", round(model.inner.incumbent.objective_value, digits=4))
    end 

    if length(model.inner.redundant_constraint_ids.linear_constraint_ids) > 0 
        @printf("%-60s : %5d\n", "Number of redundant linear constraints", length(model.inner.redundant_constraint_ids.linear_constraint_ids))
    end 

    if length(model.inner.redundant_constraint_ids.quadratic_constraint_ids) > 0 
        @printf("%-60s : %5d\n", "Number of redundant quadratic constraints", length(model.inner.redundant_constraint_ids.quadratic_constraint_ids))
    end 

    if length(model.inner.redundant_constraint_ids.soc_constraint_ids) > 0 
        @printf("%-60s : %5d\n", "Number of redundant SOC constraints", length(model.inner.redundant_constraint_ids.soc_constraint_ids))
    end 

    if length(model.inner.redundant_constraint_ids.rsoc_constraint_ids) > 0 
        @printf("%-60s : %5d\n", "Number of redundant RSOC constraints", length(model.inner.redundant_constraint_ids.rsoc_constraint_ids))
    end 

    if length(model.inner.redundant_constraint_ids.nl_constraint_ids) > 0 
        @printf("%-60s : %5d\n", "Number of redundant nonlinear constraints", length(model.inner.redundant_constraint_ids.nl_constraint_ids))
    end 


    return

end