"""
	Tree traversal for bi-linear term detection
	Input: expression subtree | terms::Array | buffer::Array
	Output: the non-linear terms and an empty buffer
	Note: buffer delimits the expression on :(+) and :(-)

	Algorithm: DFS tree traversal to collect bi-linear terms.

	NOTE: This function should be modifined when working on with multi-linears
"""
function traverse_expr(expr, terms::Array=[], buffer::Array=[]; kwargs...)


	if isa(expr, Float64) || isa(expr, Int64)
		return terms, buffer
	elseif expr == :+ || expr == :-
		if !(buffer in terms) && !isempty(buffer)
			push!(terms, buffer)
		end
		buffer = []
		return terms, buffer
	elseif expr == :*
		return terms, buffer
	elseif expr == :^
		return terms, buffer
	elseif expr == :/
		error("Unspported operator $(expr)")
	elseif expr.head == :ref
		push!(buffer, expr)
		return terms, buffer
	end

	for i in 1:length(expr.args)
		terms, buffer = traverse_expr(expr.args[i], terms, buffer)

		if expr.args[1] == :+ || expr.args[1] == :-
			if !(buffer in terms) && !isempty(buffer)
				push!(terms, buffer)
			end
			buffer = []
		end
	end

	# handling special operators :(^)
	if expr.args[1] == :^
		push!(buffer, buffer[1])
		if !(buffer in terms) && !isempty(buffer)
			push!(terms, buffer)
		end
		buffer = []
	end

	return terms, buffer
end

"""
	Populate the dictionary that contains the information of all the non-linear terms
"""
function populate_dict_nonlinear_info(m::PODNonlinearModel; kwargs...)

	# options is not used anywhere
	# options = Dict(kwargs)

	m.lifted_obj_expr_mip = deepcopy(m.obj_expr_orig)
	m.lifted_constr_expr_mip = []
 	for i in 1:m.num_constr_orig
		push!(m.lifted_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
	end

	m.dict_nonlinear_info[:xdim] = m.num_var_orig
	terms = []
	analyze_expr(m.lifted_obj_expr_mip, terms, m.dict_nonlinear_info)

	for i in 1:m.num_constr_orig
		analyze_expr(m.lifted_constr_expr_mip[i].args[2], terms, m.dict_nonlinear_info)
	end

	return m
end


"""
	Dereferencing function: splice and dereference
"""
function expr_dereferencing()

end

"""
	Analyze a single expression and create lifted terms and in-place lifted expression

	Input: expression | existing terms already parsed | dict_nonlinear_info
	Output: Updated terms array and dict_nonlinear_info (not returned, updated data structures)
	Algorithm:
		1) traverse the expression for term detection using traverse_expr
		2) lift the original expression in place

"""
function analyze_expr(expr, terms::Array=[], dict_nonlinear_info::Dict=Dict(); kwargs...)

	# options not being used anywhere
	# options = Dict(kwargs)
	terms, buffer = traverse_expr(expr, terms)
	add_info_using_terms(terms, dict_nonlinear_info)

	return
end


"""
	Adding the lifted variable map using the non-linear terms

	Input: array of non-linear terms | dict_nonlinear_info
	Output: updated dict_nonlinear_info (no explicit return)

	Mapping : term => Dict(:lifted_var_ref::Expr, :key::Expr, :lifted_constr_ref::Expr)

"""
function add_info_using_terms(terms::Array=[], dict_nonlinear_info=Dict(); kwargs...)

	@assert haskey(dict_nonlinear_info, :xdim)
	yidx = dict_nonlinear_info[:xdim] + length(dict_nonlinear_info)
	# Logic only works for bi-linear terms
	for term in terms
		if length(term) > 1 && !(term in keys(dict_nonlinear_info)) && !(reverse(term) in keys(dict_nonlinear_info)) # Check for duplicates
			lifted_var_ref = Expr(:ref, :x, yidx)
            lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, term[1], term[2]))
			dict_nonlinear_info[term] = Dict(:lifted_var_ref => lifted_var_ref, :ref => term, :lifted_constr_ref => lifted_constr_ref)
			yidx += 1
		end
	end

	return
end

# from here on this code was written by sitew and has not been verified.

"""
	Replace expression variable with lifted variable in-place
	Input: expression::Expr | dict_nonlinear_info::Dict
	Output: recursive in-place operation
"""
function expr_lift(expr, dict_nonlinear_info::Dict; kwargs...)

	@assert expr_mark(deepcopy(expr)) <= 2 # Not focusing on multi-linear terms
	expr_flatten(expr)	# Preproces the expression
	if expr.args[1] in [:+,:-]  || !isempty([i for i in 1:length(expr.args) if (isa(expr.args[i], Expr) && expr.args[i].head == :call)])
		for i in 2:length(expr.args) # Reserved 1 with operator
			if isa(expr.args[i], Expr) && expr.args[i].head == :call
				rep = expr_lift(expr.args[i], dict_nonlinear_info) # Continue with sub-tree
				expr.args[i] = rep
			elseif isa(expr.args[i], Float64) || expr.args[i].head == :ref
				expr.args[i] = expr.args[i]
			end
		end

	elseif expr.args[1] == :* # with assumption :: bilinear, limited () usage
		idxs = [i for i in 1:length(expr.args) if (isa(expr.args[i], Expr) && expr.args[i].head == :ref)]
		refs = expr.args[idxs]
		if length(refs) == length(expr.args) - 1  # (*)->[(x),(x)] ==> (*)->[(1),(y)]
			push!(expr.args, 1.0)
		end
		if haskey(dict_nonlinear_info, refs)
			deleteat!(expr.args, idxs)
			push!(expr.args, dict_nonlinear_info[refs][:lifted_var_ref]) # In-place Lift
		elseif haskey(dict_nonlinear_info, reverse(refs)) # No duplicates
			deleteat!(expr.args, idxs)
			push!(expr.args, dict_nonlinear_info[reverse(refs)][:lifted_var_ref]) # Place Lift
		end

	elseif expr.args[1] == :^  # with assumption :: square
		refs = [v for v in expr.args if (isa(v, Expr) && v.head == :ref)]
		power = [v for v in expr.args if (isa(v, Float64) || isa(v, Int64))]
		@assert length(refs) == 1  # Error handling
		@assert length(power) == 1
		@assert power[1] == 2
		refs = [refs;refs]
		@assert haskey(dict_nonlinear_info, refs)
		expr = dict_nonlinear_info[refs][:lifted_var_ref] # Lift
	else
		error("Unspported operator for variable lifting $(expr.args[1])")
	end

	return expr
end

"""
	Recursivly pre-process expression by treating the coefficients
	TODO: this function requires a lot of refining.
	Most issues can be caused by this function.
"""
function expr_flatten(expr, level=0; kwargs...)

	if level > 0  # No trivial constraint is allowed "3>5"
		flat = expr_arrangeargs(expr.args)
		if isa(flat, Float64)
			return flat
		else
			expr.args = flat
		end
	end

	for i in 2:length(expr.args)
		if isa(expr.args[i], Expr) && expr.args[i].head == :call
			expr.args[i] = expr_flatten(expr.args[i], level+1)
		end
	end
	if level > 0 #Root level process, no additional operators
		expr.args = expr_arrangeargs(expr.args)
	end

	return expr
end

"""
	Re-arrange childrens by their type.
	Only consider 3 types: coefficients(Int64, Float64), :ref, :calls
	Sequence arrangement is dependent on operator
"""
function expr_arrangeargs(args::Array; kwargs...)

	# A mapping used for operator treatments
	reverter = Dict(true=>:-, false=>:+)
	operator_coefficient_base = Dict{Symbol, Float64}(
		:+ => 0.0,
		:- => 0.0,
		:/ => 1.0,
		:* => 1.0,
		:^ => 1.0
	)

	# Treatment => Flatten pure number sub-tree
	allFloats = [i for i in 2:length(args) if isa(args[i], Float64)]
	if length(allFloats) == length(args)-1
		val = args[2]
		for i in 3:length(args)
			val = eval(args[1])(val, args[i])
		end
		return val
	end

	if args[1] in [:^]
		return args
	elseif args[1] in [:/]
		error("Unspported operator $(args[1])")
	elseif args[1] in [:*, :+, :-]

		val = operator_coefficient_base[args[1]] # initialize
		exprs = []
		refs = []

		valinit = false  # what is first children : needs consideration with :(-)
		refinit = false
		callinit = false

		for i in 2:length(args)
			if isa(args[i],Float64)
				valinit += (i==2)
				val = eval(args[1])(val, eval(reverter[(i==2) && (args[1]==:-)])(args[i])) # Revert the first sigh if (-)
			elseif typeof(args[i]) == Expr
				if args[i].head == :ref
					refinit += (i==2)
					push!(refs, args[i])
				elseif args[i].head == :call
					callinit += (i==2)
					push!(exprs, args[i])
				else
					error("Unkown expression head")
				end
			else
				error("Unkown argument type")
			end
		end

		# Treatment for 0.0 coefficients
		if val == operator_coefficient_base[args[1]] && val == 0.0
			val = []
		end

		# Re-arrange children :: without actually doing anything on them
		if Bool(valinit)
			return [args[1];val;refs;exprs]
		elseif Bool(refinit)
			return [args[1];refs;val;exprs]
		elseif Bool(callinit)
			return [args[1];exprs;refs;val]
		else
			error("Unexpected condition encountered...")
		end

	else
		error("Unkown operator $(args[1])")
	end

end


"""
	Mark expression with "dimensions" of the variable. Augment the tree.
"""
function expr_mark(expr; kwargs...)

	p = 0
	for i in 2:length(expr.args)
		if isa(expr.args[i], Expr) && expr.args[i].head == :ref
			if expr.args[1] in [:+, :-]
				p = max(p, 1)
				expr.args[i].typ = 1
			elseif expr.args[1] in [:*]
				p = +(p, 1)
				expr.args[i].typ = 1
			elseif expr.args[1] in [:^] # Consider leaf when encounter :^
				p = expr_dim_enquiry(expr)
			elseif expr.args[1] in [:/]
				error("Does not support :/ and :^ yet.")
			end
		elseif isa(expr.args[i], Expr) && expr.args[i].head == :call
			if expr.args[1] in [:+,:-]
				innerp = expr_mark(expr.args[i])
				p = max(p, innerp)
				expr.args[i].typ = innerp
			elseif expr.args[1] in [:*]
				innerp = expr_mark(expr.args[i])
				p = +(p, innerp)
				expr.args[i].typ = innerp
			elseif expr.args[1] in [:/, :^]
				error("Does not support :/ and :^ yet.")
			end
		elseif isa(expr.args[i],Float64)
			p = p
		else
			error("Unexpected expression node encouter $(expr.args[i])")
		end
	end

	return p
end

"""
	A special check when encouter :^. Cleaner function
"""
function expr_dim_enquiry(expr)

	@assert expr.args[1] == :^
	@assert length(expr.args) == 3 # Don't take cray :call with :^
	if isa(expr.args[2], Float64)
		return expr.args[2]
	else
		return expr.args[3]
	end
end

"""
	Check if a sub-tree is linear or not
"""
function expr_islinear(expr)
	if expr_mark(expr) > 1
		return False
	else
		return True
	end
end
