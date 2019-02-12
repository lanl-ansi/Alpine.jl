
mutable struct Incumbent 
    variable_value                  ::Vector{Float64}
    objective_value                 ::Float64 
    best_bound                      ::Float64 
    status                          ::Symbol 
end

mutable struct Status 
    local_solve_status              ::Symbol 
    mip_solve_status                ::Symbol 
    alpine_status                   ::Symbol 
    bt_status                       ::Symbol 
end 


mutable struct AlpineProblem 
    # variable and constraint count
    num_variables                   ::Int 
    num_linear_le_constraints       ::Int 
    num_linear_ge_constraints       ::Int 
    num_linear_eq_constraints       ::Int  
    num_quadratic_le_constraints    ::Int  
    num_quadratic_ge_constraints    ::Int  
    num_quadratic_eq_constraints    ::Int  
    num_soc_constraints             ::Int  
    num_rsoc_constraints            ::Int  
    num_nlp_constraints             ::Int  
    
    # JuMP models 
    mip                             ::JuMP.model 
    continuous_relaxation           ::JuMP.model 

    # Variable bounds information 
    lower_original                  ::Vector{Float64}
    upper_original                  ::Vector{Float64}
    lower_tightened                 ::Vector{Float64}
    upper_tightened                 ::Vector{Float64} 

    # Nonlinear information 
    nl_terms                        ::Dict{Any,Any}

    # Incumbent information 
    incumbent                       ::Incumbent

    # Algorithm status
    status                          ::Status

    # Algorithm progress 
    iteration                       ::Int

end 

