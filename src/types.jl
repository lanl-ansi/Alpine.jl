
mutable struct Incumbent 
    variable_value                  ::Vector{Float64}
    objective_value                 ::Float64 
    best_bound                      ::Float64 
    status                          ::Symbol 
end

Incumbent() = Incumbent(Vector{Float64}(), NaN, NaN, :NaN)

mutable struct Status 
    local_solve_status              ::Symbol 
    mip_solve_status                ::Symbol 
    alpine_status                   ::Symbol 
    bt_status                       ::Symbol 
end 

Status() = Status(:NaN, :NaN, :NaN, :NaN)

mutable struct SolverOptions 
    nlp_optimizer                   ::Union{Nothing, MOI.AbstractOptimizer}
    minlp_optimizer                 ::Union{Nothing, MOI.AbstractOptimizer}
    mip_optimizer                   ::Union{Nothing, MOI.AbstractOptimizer}
    
    log_level                       ::Int 
    time_limit                      ::Float64 
    max_iter                        ::Int 
    rel_gap                         ::Float64 
    abs_gap                         ::Float64 
    opt_tol                         ::Float64 
    
    disc_partition_method           ::Symbol 
    disc_ratio                      ::Float64 
    disc_uniform_rate               ::Float64 
    disc_divert_chunks              ::Int 
    disc_abs_width_tol              ::Float64 
    disc_rel_width_tol              ::Float64 
    disc_consecutive_forbid         ::Int 

    bt                              ::Bool 
    presolve_time_limit             ::Float64 
    bt_max_iter                     ::Int 
    bt_width_tol                    ::Float64
    bt_improvement_tol              ::Float64 
    bt_precision                    ::Int 
    bt_algo                         ::Symbol 
    bt_relax                        ::Bool 
    bt_mip_time_limit               ::Float64 
    bp                              ::Bool 

    is_convex                       ::Bool
end 

mutable struct AlpineProblem 
    # variable and constraint count
    num_variables                   ::Int 
    num_constraints                 ::Int
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
    mip                             ::Union{Nothing, JuMP.Model}
    continuous_relaxation           ::Union{Nothing, JuMP.Model} 

    # Variable bounds information 
    lower_original                  ::Union{Nothing, Vector{Float64}}
    upper_original                  ::Union{Nothing, Vector{Float64}}
    lower_tightened                 ::Union{Nothing, Vector{Float64}}
    upper_tightened                 ::Union{Nothing, Vector{Float64}}

    # Nonlinear information 
    is_objective_linear             ::Union{Nothing, Bool} 
    is_objective_quadratic          ::Union{Nothing, Bool} 
    is_objective_nl                 ::Union{Nothing, Bool} 
    objective_expr                  ::Union{Nothing, Expr} 
    nl_constraint_expr              ::Union{Nothing, Vector{Expr}}
    nl_terms                        ::Union{Nothing, Dict{Expr, Any}}
    constraints_with_nl_terms       ::Union{Nothing, Vector{Int}}
    lifted_constraints              ::Union{Nothing, Vector{JuMP.ConstraintRef}}
    lifted_var_info                 ::Union{Nothing, Dict{Expr, Any}}

    # Incumbent information 
    incumbent                       ::Union{Nothing, Incumbent}

    # Algorithm status
    status                          ::Union{Nothing, Status}

    # Algorithm progress 
    iteration                       ::Int

end 

AlpineProblem() = AlpineProblem(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    nothing, nothing, 
    nothing, nothing, nothing, nothing, 
    nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
    nothing, 
    nothing, 
    0
)
     