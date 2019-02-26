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
    
    num_linear_constraints = ap.num_linear_eq_constraints + 
        ap.num_linear_ge_constraints + ap.num_linear_le_constraints 
    num_quadratic_constraints = ap.num_quadratic_eq_constraints +
        ap.num_quadratic_ge_constraints + ap.num_quadratic_le_constraints 
    num_soc_constraints = ap.num_soc_constraints 
    num_rsoc_constraints = ap.num_rsoc_constraints 

    @printf("%-40s : %5d\n", "Number of constraints", ap.num_constraints)
    
    if num_linear_constraints != 0
        @printf("%-40s : %5d\n", "Number of linear constraints", num_linear_constraints)
    end 

    if num_quadratic_constraints != 0 
        @printf("%-40s : %5d\n", "Number of quadratic constraints", num_quadratic_constraints)
    end 

    if num_soc_constraints != 0
        @printf("%-40s : %5d\n", "Number of SOC constraints", num_soc_constraints)
    end 

    if num_rsoc_constraints != 0 
        @printf("%-40s : %5d\n", "Number of RSOC constraints", num_rsoc_constraints)
    end

    if ap.num_nlp_constraints != 0
        @printf("%-40s : %5d\n", "Number of nonlinear constraints", ap.num_nlp_constraints)
    end 

    println()
end