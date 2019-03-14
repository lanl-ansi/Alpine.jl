"""
MOI variables
"""
    
"""
Getters for variables 
"""
MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)

function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i in 1:length(model.variable_info)]
end 

"""
``MOI.add_variable()`` overloads and safety functions
"""
function MOI.add_variable(model::Optimizer)
    push!(model.variable_info, VariableInfo())
    return MOI.VariableIndex(length(model.variable_info))
end

function MOI.add_variables(model::Optimizer, n::Int)
    return [MOI.add_variable(model) for i in 1:n]
end

function check_inbounds(model::Optimizer, vi::VI)
	num_variables = length(model.variable_info)
	if !(1 <= vi.value <= num_variables)
	    @error "Invalid variable index $vi. ($num_variables variables in the model.)"
	end
	return
end
	
check_inbounds(model::Optimizer, var::SVF) = check_inbounds(model, var.variable)
	
function check_inbounds(model::Optimizer, aff::SAF)
	for term in aff.terms
	    check_inbounds(model, term.variable_index)
	end
	return
end
	
function check_inbounds(model::Optimizer, quad::SQF)
	for term in quad.affine_terms
	    check_inbounds(model, term.variable_index)
	end
	for term in quad.quadratic_terms
	    check_inbounds(model, term.variable_index_1)
	    check_inbounds(model, term.variable_index_2)
	end
	return
end
	
has_upper_bound(model::Optimizer, vi::VI) = 
    model.variable_info[vi.value].has_upper_bound
	
has_lower_bound(model::Optimizer, vi::VI) = 
    model.variable_info[vi.value].has_lower_bound
	
is_fixed(model::Optimizer, vi::VI) = 
    model.variable_info[vi.value].is_fixed

is_binary(model::Optimizer, vi::VI) = 
    model.variable_info[vi.value].is_binary

""" 
Primal-start support for Alpine
"""
MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{VI}) = true

function MOI.set(model::Optimizer, ::MOI.VariablePrimalStart,
    vi::VI, value::Union{Real, Nothing})
    check_inbounds(model, vi)
    model.variable_info[vi.value].start = value
    return
end

"""
MOI.VariableName attribute support 
"""
MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{VI}) = true

function MOI.set(model::Optimizer, ::MOI.VariableName, vi::VI, name::String)
	model.variable_info[vi.value].name = name
	return
end

function MOI.get(model::Optimizer, ::MOI.VariableName, vi::VI)
	return model.variable_info[vi.value].name
end

"""
Populating Variable bounds 
"""
function MOI.add_constraint(model::Optimizer, v::SVF, lt::MOI.LessThan{Float64})
    vi = v.variable
    check_inbounds(model, vi)
    if isnan(lt.upper)
        @error "Invalid upper bound value $(lt.upper)."
    end
    if has_upper_bound(model, vi)
        @error "Upper bound on variable $vi already exists."
    end
    if is_fixed(model, vi)
        @error "Variable $vi is fixed. Cannot also set upper bound."
    end
    model.variable_info[vi.value].upper_bound = lt.upper
    model.variable_info[vi.value].has_upper_bound = true
    return MOI.ConstraintIndex{SVF, MOI.LessThan{Float64}}(vi.value)
end

function MOI.add_constraint(model::Optimizer, v::SVF, gt::MOI.GreaterThan{Float64})
    vi = v.variable
    check_inbounds(model, vi)
    if isnan(gt.lower)
        @error "Invalid lower bound value $(gt.lower)."
    end
    if has_lower_bound(model, vi)
        @error "Lower bound on variable $vi already exists."
    end
    if is_fixed(model, vi)
        @error "Variable $vi is fixed. Cannot also set lower bound."
    end
    model.variable_info[vi.value].lower_bound = gt.lower
    model.variable_info[vi.value].has_lower_bound = true
    return MOI.ConstraintIndex{SVF, MOI.GreaterThan{Float64}}(vi.value)
end

function MOI.add_constraint(model::Optimizer, v::SVF, eq::MOI.EqualTo{Float64})
    vi = v.variable
    check_inbounds(model, vi)
    if isnan(eq.value)
        @error "Invalid fixed value $(gt.lower)."
    end
    if has_lower_bound(model, vi)
        @error "Variable $vi has a lower bound. Cannot be fixed."
    end
    if has_upper_bound(model, vi)
        @error "Variable $vi has an upper bound. Cannot be fixed."
    end
    if is_fixed(model, vi)
        @error "Variable $vi is already fixed."
    end
    model.variable_info[vi.value].lower_bound = eq.value
    model.variable_info[vi.value].upper_bound = eq.value
    model.variable_info[vi.value].is_fixed = true
    model.variable_info[vi.value].is_bounded = true
    return MOI.ConstraintIndex{SVF, MOI.EqualTo{Float64}}(vi.value)
end
