"""
    matrix_from_quadratic_terms(terms::Vector{MOI.ScalarQuadraticTerm})

Creates a sparse matrix from quadratic terms and returns a Q and index_to_variable_map
"""
function matrix_from_quadratic_terms(terms::Vector{SQT})
    variable_to_index_map = Dict{VI, Int}()
    index_to_variable_map = Dict{Int, VI}()
    n = 0
    for term in terms
        for variable in (term.variable_index_1, term.variable_index_2)
            if !(variable in keys(variable_to_index_map))
                n += 1
                variable_to_index_map[variable] = n
                index_to_variable_map[n] = variable
            end
        end
    end

    I = Int[]
    J = Int[]
    V = Float64[]
    
    for term in terms
        i = variable_to_index_map[term.variable_index_1]
        j = variable_to_index_map[term.variable_index_2]
        push!(I, i)
        push!(J, j)
        push!(V, term.coefficient)
        if i != j
            push!(I, j)
            push!(J, i)
            push!(V, term.coefficient)
        end
    end

    # Duplicate terms are summed together in `sparse`
    Q = sparse(I, J, V, n, n)
    
    return Q, index_to_variable_map
end

"""
    matrix_from_nl_function(nl_function::NLFunction)

Creates a sparse matrix from the quadratic part of an nonlinear function data object and returns a Q and index_to_variable_map
"""

function matrix_from_nl_function(nl_function::NLFunction)
    quadratic_part = nl_function.quadratic_part 
    bilinear_part = nl_function.bilinear_part 
    if isa(quadratic_part, Nothing) && isa(bilinear_part, Nothing)
        return nothing, nothing 
    end
    variable_to_index_map = Dict{VI, Int}()
    index_to_variable_map = Dict{Int, VI}()
    n = 0
    if ~isa(quadratic_part, Nothing)
        for quad_term in quadratic_part 
            coeff, expr = quad_term.expression 
            variable = expr.args[2].args[2]
            if !(variable in keys(variable_to_index_map))
                n += 1
                variable_to_index_map[variable] = n
                index_to_variable_map[n] = variable
            end
        end 
    end 

    if ~isa(bilinear_part, Nothing)
        for bilinear_term in bilinear_part 
            coeff, expr = bilinear_term.expression 
            variable_1 = expr.args[2].args[2]
            variable_2 = expr.args[3].args[2]
            for variable in (variable_1, variable_2)
                if !(variable in keys(variable_to_index_map))
                    n += 1
                    variable_to_index_map[variable] = n
                    index_to_variable_map[n] = variable
                end
            end
        end
    end 

    I = Int[]
    J = Int[]
    V = Float64[]

    if ~isa(quadratic_part, Nothing)
        for quad_term in quadratic_part 
            coeff, expr = quad_term.expression 
            variable = expr.args[2].args[2]
            i = variable_to_index_map[variable]
            push!(I, i)
            push!(J, i)
            push!(V, coeff * 2.0)
        end 
    end 

    if ~isa(bilinear_part, Nothing)
        for bilinear_term in bilinear_part 
            coeff, expr = bilinear_term.expression 
            variable_1 = expr.args[2].args[2]
            variable_2 = expr.args[3].args[2]
            i = variable_to_index_map[variable_1]
            j = variable_to_index_map[variable_2]
            push!(I, i)
            push!(J, j)
            push!(V, coeff)
            push!(I, j)
            push!(J, i)
            push!(V, coeff)
        end
    end 

    # Duplicate terms are summed together in `sparse`
    Q = sparse(I, J, V, n, n)
    
    return Q, index_to_variable_map
end 

"""
    get_type_dict(obj::Dict)::Dict{Symbol,Type}() 

This function returns the types of every field in a given dictionary 
"""
function get_type_dict(obj)
    T = typeof(obj)
    type_dict = Dict{Symbol,Type}()
    for (name, typ) in zip(fieldnames(T), T.types)
        type_dict[name] = typ
    end
    return type_dict
end

"""
    dereference_expr!(expr, m::JuMP.model)

This function replaces the MOI variables in an expression tree using JuMP VariableRef 
"""
function dereference_expr!(expr, m)
    for i in 2:length(expr.args)
        if isa(expr.args[i], Union{Float64,Int64})
            continue
        elseif expr.args[i].head == :ref
            @assert isa(expr.args[i].args[2], MOI.VariableIndex)
            expr.args[i] = JuMP.VariableRef(m, expr.args[i].args[2])
        elseif expr.args[i].head == :call
            dereference_expr!(expr.args[i], m)
        else
            error(LOGGER, "dereference_expr :: Unexpected term in expression tree")
        end
    end

    return
end

"""
    state_is_optimal(state::MOI.TerminationStatusCode; allow_almost=false)

Returns true if either status is optimal or locally solved. If allow_almost=true then check for `ALMOST_LOCALLY_SOLVED`
"""
function status_is_optimal(status::MOI.TerminationStatusCode; allow_almost=false)
    return status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || (allow_almost && status == MOI.ALMOST_LOCALLY_SOLVED)
end

"""
    clear_optimizers!(model::MOI.AbstractOptimizer)

Drops the optimizers from the models so that resolves can be peformed without re-initializing Alpine model 
"""
function clear_optimizers!(model::MOI.AbstractOptimizer)
    if ~isa(model.inner.mip, Nothing)
        if JuMP.backend(model.inner.mip).state == MOIU.ATTACHED_OPTIMIZER
            MOIU.reset_optimizer(JuMP.backend(model.inner.mip))
        end
    end 
    if ~isa(model.inner.continuous_relaxation, Nothing)
        if JuMP.backend(model.inner.continuous_relaxation).state == MOIU.ATTACHED_OPTIMIZER
            MOIU.reset_optimizer(JuMP.backend(model.inner.continuous_relaxation))
        end
    end 
    if ~isa(model.inner.obbt_model, Nothing)
        if JuMP.backend(model.inner.obbt_model).state == MOIU.ATTACHED_OPTIMIZER
            MOIU.reset_optimizer(JuMP.backend(model.inner.obbt_model))
        end
    end

    return
end 