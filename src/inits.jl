"""
init_ap_data!() for Alpine 
""" 
function init_ap_data!(model::MOI.AbstractOptimizer)
    model.inner = AlpineProblem()
    model.inner.num_variables = length(model.variable_info)
    model.inner.num_linear_le_constraints = length(model.linear_le_constraints)
    model.inner.num_linear_ge_constraints = length(model.linear_ge_constraints)
    model.inner.num_linear_eq_constraints = length(model.linear_eq_constraints)
    model.inner.num_quadratic_le_constraints = length(model.quadratic_le_constraints)
    model.inner.num_quadratic_ge_constraints = length(model.quadratic_ge_constraints)
    model.inner.num_quadratic_eq_constraints = length(model.quadratic_eq_constraints)
    model.inner.num_soc_constraints = length(model.soc_constraints)
    model.inner.num_rsoc_constraints = length(model.rsoc_constraints) 

    model.inner.is_constraint_convex = Dict{Symbol, Vector{Bool}}()

    if model.inner.num_quadratic_le_constraints > 0 
        model.inner.is_constraint_convex[:quadratic_le] = 
            [false for i in 1:model.inner.num_quadratic_le_constraints]
    end 

    if model.inner.num_quadratic_ge_constraints > 0 
        model.inner.is_constraint_convex[:quadratic_ge] = 
            [false for i in 1:model.inner.num_quadratic_ge_constraints]
    end

    model.inner.lower_original = Vector{Float64}()
    model.inner.upper_original = Vector{Float64}()

    for vi in model.variable_info
        if !vi.has_lower_bound || !vi.has_upper_bound 
            error(LOGGER, "Alpine.jl requires every variable in the problem to be bounded; ensure that finite bounds are provided for every variable using JuMP.set_lower_bound() and JuMP.set_upper_bound()")
        end 

        push!(model.inner.lower_original, vi.lower_bound)
        push!(model.inner.upper_original, vi.upper_bound)
    end 

    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)
        evaluator = model.nlp_data.evaluator
        MOI.initialize(evaluator, [:ExprGraph])
        num_nlp_constraints = length(model.nlp_data.constraint_bounds)
        model.inner.num_nlp_constraints = num_nlp_constraints
        if num_nlp_constraints > 0 
            model.inner.nl_constraint_expr = Vector{Expr}()
        end
        
        if model.nlp_data.has_objective 
            model.inner.is_objective_nl = true 
            model.inner.objective_expr = MOI.objective_expr(evaluator)
            model.inner.is_objective_linear = false 
            model.inner.is_objective_quadratic = false 
        elseif isa(model.objective, SQF)
            model.inner.is_objective_nl = false
            model.inner.is_objective_linear = false 
            model.inner.is_objective_quadratic = true
        elseif isa(model.objective, Union{SAF, SVF})
            model.inner.is_objective_nl = false
            model.inner.is_objective_linear = true
            model.inner.is_objective_quadratic = false
            model.is_objective_convex = true
        end

        for i in 1:num_nlp_constraints
            constraint_expr = MOI.constraint_expr(evaluator, i)
            push!(model.inner.nl_constraint_expr, constraint_expr)
        end

    else 
        info(LOGGER, "no explicit NLP constraints or objective provided using @NLconstraint or @NLobjective macros")
        model.inner.num_nlp_constraints = 0
        model.inner.is_objective_nl = false
        if isa(model.objective, SQF)
            model.inner.is_objective_quadratic = true 
            model.inner.is_objective_linear = false 
        elseif isa(model.objective, Union{SAF, SVF})
            model.inner.is_objective_linear = true 
            model.inner.is_objective_quadratic = false 
            model.is_objective_convex = true
        end
    end 

    model.inner.num_constraints = 
        model.inner.num_linear_le_constraints + 
        model.inner.num_linear_ge_constraints +
        model.inner.num_linear_eq_constraints +
        model.inner.num_quadratic_le_constraints +
        model.inner.num_quadratic_ge_constraints +
        model.inner.num_quadratic_eq_constraints +
        model.inner.num_soc_constraints +
        model.inner.num_rsoc_constraints + 
        model.inner.num_nlp_constraints
    
    if model.solver_options.log_level != 0
        print_var_con_summary(model.inner)
    end
    
end