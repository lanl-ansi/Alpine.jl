"""
init_ap_data!() for Alpine 
""" 
function init_ap_data!(model::MOI.AbstractOptimizer)
    model.inner = AlpineProblem()
    # populate variable/constraint statistics
    model.inner.num_variables = length(model.variable_info)
    model.inner.num_linear_constraints = length(model.linear_constraints)
    model.inner.num_quadratic_constraints = length(model.quadratic_constraints)
    model.inner.num_soc_constraints = length(model.soc_constraints)
    model.inner.num_rsoc_constraints = length(model.rsoc_constraints) 

    # initialize constraint bound vectors
    model.inner.constraint_bound_info = Vector{Interval{Float64}}()
    model.inner.objective_bound_info = -Inf..Inf

    # populate linear constraint bounds
    for (func, set) in model.linear_constraints 
        if isa(set, MOI.LessThan{Float64}) 
            push!(model.inner.constraint_bound_info, -Inf..set.upper)
        elseif isa(set, MOI.GreaterThan{Float64})
            push!(model.inner.constraint_bound_info, set.lower..Inf)
        else 
            push!(model.inner.constraint_bound_info, set.value..set.value)
        end 
    end 

    # populate quadratic constraint bounds
    for (func, set) in model.quadratic_constraints 
        if isa(set, MOI.LessThan{Float64}) 
            push!(model.inner.constraint_bound_info, -Inf..set.upper)
        elseif isa(set, MOI.GreaterThan{Float64})
            push!(model.inner.constraint_bound_info, set.lower..Inf)
        else 
            push!(model.inner.constraint_bound_info, set.value..set.value)
        end 
    end 

    # handle SOCs and RSOCs later with additional book-keeping in AlpineProblem struct
    for constraint in model.soc_constraints 
        push!(model.inner.constraint_bound_info, -Inf..Inf)
    end 

    for constraint in model.rsoc_constraints 
        push!(model.inner.constraint_bound_info, -Inf..Inf)
    end

    # initialize convexity vectors and booleans
    model.inner.quadratic_constraint_convexity = Symbol[]
    model.inner.quadratic_function_convexity = Symbol[]

    model.inner.is_objective_linear = false 
    model.inner.is_objective_nl = false 
    model.inner.is_objective_quadratic = false

    if model.inner.num_quadratic_constraints > 0 
        model.inner.quadratic_constraint_convexity = 
            [:undet for i in 1:model.inner.num_quadratic_constraints]
        model.inner.quadratic_function_convexity = 
            [:undet for i in 1:model.inner.num_quadratic_constraints]
    end 
    
    # initialize variable bound vectors
    model.inner.variable_bound_original = Vector{Interval{Float64}}()
    model.inner.variable_bound_tightened = Vector{Interval{Float64}}()

    # populate number of binary and integer variables 
    binary_variables = filter(i -> i == true, info_array_of_variables(model.variable_info, :is_binary))
    model.inner.num_binary_variables = length(binary_variables)


    # populate variable bound vectors
    for vi in model.variable_info
        lb = -Inf 
        ub = Inf
        (vi.has_lower_bound) && (lb = vi.lower_bound)
        (vi.has_upper_bound) && (ub = vi.upper_bound)

        push!(model.inner.variable_bound_original, lb..ub)
        push!(model.inner.variable_bound_tightened, lb..ub)
    end 

    # if NLP populate NLP part 
    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)
        evaluator = model.nlp_data.evaluator
        MOI.initialize(evaluator, [:ExprGraph])
        num_nlp_constraints = length(model.nlp_data.constraint_bounds)
        model.inner.num_nlp_constraints = num_nlp_constraints
        if num_nlp_constraints > 0 
            model.inner.nl_constraint_expression = Vector{Expr}()
        end

        for i in 1:num_nlp_constraints 
            lb = model.nlp_data.constraint_bounds[i].lower
            ub = model.nlp_data.constraint_bounds[i].upper
            push!(model.inner.constraint_bound_info, lb..ub)
        end 
        
        if model.nlp_data.has_objective 
            model.inner.is_objective_nl = true 
            model.inner.objective_expression = MOI.objective_expr(evaluator)
            model.inner.is_objective_linear = false 
            model.inner.is_objective_quadratic = false 
            model.inner.objective_convexity = :undet
        elseif isa(model.objective, SQF)
            model.inner.is_objective_nl = false
            model.inner.is_objective_linear = false 
            model.inner.is_objective_quadratic = true
            model.inner.objective_convexity = :undet
        elseif isa(model.objective, Union{SAF, SVF})
            model.inner.is_objective_nl = false
            model.inner.is_objective_linear = true
            model.inner.is_objective_quadratic = false
            model.inner.objective_convexity = :convex
        end

        for i in 1:num_nlp_constraints
            constraint_expr = MOI.constraint_expr(evaluator, i)
            push!(model.inner.nl_constraint_expression, constraint_expr)
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
            model.objective_convexity = true
        end
    end 

    # compute total number of constraints 
    model.inner.num_constraints = 
        model.inner.num_linear_constraints + 
        model.inner.num_quadratic_constraints +
        model.inner.num_soc_constraints +
        model.inner.num_rsoc_constraints + 
        model.inner.num_nlp_constraints
    
    # print variable/constraint summary
    if model.solver_options.log_level != 0
        print_var_con_summary(model.inner)
    end 

    # create NLFunction for quadratic constraints 
    create_quadratic_nl_functions!(model)

    # create NLFunction for nonlinear constraints 
    clean_nl_expressions!(model)
    nl_function = NLFunction[]
    
    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)
        for i in 1:model.inner.num_nlp_constraints 
            expr = model.inner.nl_constraint_expression[i]
            disaggregated_expr = expr_disaggregate(expr)
            push!(nl_function, create_nl_function(disaggregated_expr))
        end
        model.inner.nl_function = nl_function

        if model.inner.is_objective_nl 
            disaggregated_expr = expr_disaggregate(model.inner.objective_expression)
            model.inner.objective_nl_function = create_nl_function(disaggregated_expr)
        end 
    end

    # create a DAG based using all the quadratic and nonlinear constraints
    create_dag!(model)

    # Compute and populate quadratic matrices for the quadratic and nonlinear constraints; quadratic/nonlinear objective 
    if model.inner.num_quadratic_constraints > 0 
        model.inner.quadratic_matrix = Vector{QuadraticMatrixInfo}() 
        for i in 1:model.inner.num_quadratic_constraints
            func = model.quadratic_constraints[i][1]
            Q, index_to_variable_map = matrix_from_quadratic_terms(func.quadratic_terms)
            push!(model.inner.quadratic_matrix, QuadraticMatrixInfo(Q, index_to_variable_map))
        end 
    end 

    if model.inner.num_nlp_constraints > 0
        model.inner.quadratic_matrix_nl = Vector{Union{QuadraticMatrixInfo, Nothing}}()
        for i in 1:model.inner.num_nlp_constraints
            func = model.inner.nl_function[i]
            Q, index_to_variable_map = matrix_from_nl_function(func)
            if isa(Q, Nothing) 
                push!(model.inner.quadratic_matrix_nl, nothing)
            else 
                push!(model.inner.quadratic_matrix_nl, QuadraticMatrixInfo(Q, index_to_variable_map))
            end
        end
    end

    if model.inner.is_objective_quadratic 
        func = model.objective
        Q, index_to_variable_map = matrix_from_quadratic_terms(func.quadratic_terms)
        model.inner.quadratic_matrix_objective = QuadraticMatrixInfo(Q, index_to_variable_map)
    elseif model.inner.is_objective_nl 
        func = model.inner.objective_nl_function 
        Q, index_to_variable_map = matrix_from_nl_function(func)
        (~isa(Q, Nothing)) && (model.inner.quadratic_matrix_objective = QuadraticMatrixInfo(Q, index_to_variable_map))
    end
    
    # build local solve model (MINLP -> continuous relaxation)


    return
end
