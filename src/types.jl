
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

    is_problem_convex               ::Bool
end 

mutable struct VariableInfo
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # Implies lower_bound == Inf
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # Implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound
    is_binary::Bool       # Implies lower_bound == 0, upper_bound == 1 and is MOI.ZeroOne
    is_bounded::Bool      # has_lower_bound == true && has_upper_bound == true 
    is_in_nl_term::Bool   # true if the variable is a part of some non-linear term 
    name::String
    start::Union{Nothing, Float64}
end

VariableInfo() = VariableInfo(-Inf, false, Inf, false, false, false, false, false, "", nothing)

function info_array_of_variables(variable_info::Vector{VariableInfo}, attr::Symbol)
    len_var_info = length(variable_info)
    type_dict = get_type_dict(variable_info[1])
    result = Array{type_dict[attr], 1}(undef, len_var_info)
    for i = 1:len_var_info
        result[i] = getfield(variable_info[i], attr)
    end
    return result
end


mutable struct TermInfo 
    lifted_variable_id              ::Int 
    convexity                       ::Symbol
    lifted_variable_info            ::VariableInfo 
    term_type                       ::Union{Nothing, Symbol}
end 

TermInfo() = TermInfo(NaN, :undet, VariableInfo(), nothing)
    
mutable struct Terms 
    quadratic_terms                 ::Union{Nothing, Dict{Expr, TermInfo}}
    bilinear_terms                  ::Union{Nothing, Dict{Expr, TermInfo}}
    multilinear_terms               ::Union{Nothing, Dict{Expr, TermInfo}}
    polynomial_terms                ::Union{Nothing, Dict{Expr, TermInfo}}
    trigonometric_terms             ::Union{Nothing, Dict{Expr, TermInfo}}
    other_terms                     ::Union{Nothing, Dict{Expr, TermInfo}}
end 

Terms() = Terms(nothing, nothing, nothing, 
    nothing, nothing, nothing)

mutable struct AlpineExpr 
    Expr                            ::Expr 
    convexity                       ::Symbol
end 

mutable struct NLFunction 
    linear_part                     ::Union{Nothing, AlpineExpr}
    quadratic_part                  ::Union{Nothing, AlpineExpr}
    polynomial_part                 ::Union{Nothing, AlpineExpr}
    other_part                      ::Union{Nothing, AlpineExpr}
end 

NLFunction() = NLFunction(nothing, nothing, nothing, nothing)


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
    nl_terms                        ::Union{Nothing, TermInfo}
    constraints_with_nl_terms       ::Union{Nothing, Vector{Int}}
    lifted_constraints              ::Union{Nothing, Vector{JuMP.ConstraintRef}}
    lifted_var_info                 ::Union{Nothing, Dict{Int, Any}}

    # convexity information 
    objective_convexity             ::Union{Nothing, Symbol} 
    quadratic_function_convexity    ::Union{Nothing, Dict{Symbol,Vector{Symbol}}}
    quadratic_constraint_convexity  ::Union{Nothing, Dict{Symbol,Vector{Symbol}}}
    nl_function_convexity           ::Union{Nothing, Dict{Symbol,Vector{Symbol}}}
    nl_constraint_convexity         ::Union{Nothing, Dict{Symbol,Vector{Symbol}}}

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
    nothing, nothing, nothing, nothing, nothing,
    nothing, 
    nothing, 
    0
)
     