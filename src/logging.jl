"""
Print line 
"""
function print_line(count)
    for i in 1:count
        print("-")
    end
    print("\n") 
end

"""
Print variable and constraint summary 
"""
function print_var_con_summary(ap::AlpineProblem)

    printstyled("\nProblem statistics\n", color=:cyan)
    @printf("%-40s : %5d\n", "Number of variables", ap.num_variables)
    
    @printf("%-40s : %5d\n", "Number of constraints", ap.num_constraints)
    
    if ap.num_linear_constraints != 0
        @printf("%-40s : %5d\n", "Number of linear constraints", ap.num_linear_constraints )
    end 

    if ap.num_quadratic_constraints != 0 
        @printf("%-40s : %5d\n", "Number of quadratic constraints", ap.num_quadratic_constraints)
    end 

    if ap.num_soc_constraints != 0
        @printf("%-40s : %5d\n", "Number of SOC constraints", ap.num_soc_constraints)
    end 

    if ap.num_rsoc_constraints != 0 
        @printf("%-40s : %5d\n", "Number of RSOC constraints", ap.num_rsoc_constraints)
    end

    if ap.num_nlp_constraints != 0
        @printf("%-40s : %5d\n", "Number of nonlinear constraints", ap.num_nlp_constraints)
    end 

    println()
end