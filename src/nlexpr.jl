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

	options = Dict(kwargs)

	m.lifted_obj_expr_mip = deepcopy(m.obj_expr_orig)
	m.lifted_constr_expr_mip = []
 	for i in 1:m.num_constr_orig
		push!(m.lifted_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
	end

	terms, m.dict_nonlinear_info = analyze_expr(m.lifted_obj_expr_mip, [], m.dict_nonlinear_info, m.num_var_orig)
	for i in 1:m.num_constr_orig
		terms, m.dict_nonlinear_info = analyze_expr(m.lifted_constr_expr_mip[i].args[2], terms, m.dict_nonlinear_info)
	end

	return m
end


"""
	Dereferencing function: splice and dereference
"""
function expr_dereferencing()

end

"""
	Main Expression Conversion Function :: for one expression
	Convert an expression to POD-like format with lifted variables in places

	Input: NLPEvaluator | Expression | Existing Terms | ML Terms Mapping
	Output: Updated Terms and Terms Mapping
	Algorithm:
		1) Expression Traversal for Term Detection
		2) Update terms and their mapping
		3) Lift original expression in-place

	Variable Convention : Operation base => MathProgBase :: x
"""
function analyze_expr(expr, term::Array=[], t2y::Dict=Dict(), colCnt::Int=0; kwargs...)

	options = Dict(kwargs)
	term, non = traverse_expr(expr, term)
	t2y = expr_t2y(term, t2y, colCnt)

	return term, t2y
end


"""
	Lift variable mapping :: initialize or continue build variable mapping

	Input: numCols::Int | lifted terms::Array | t2y mapping dictionary::Dict
	Output: updated mapping dictionary::Dict

	Naming : x[1:numCols] > [x;lifted x] > lifted x still use name x
	Mapping : [:Expr, :Expr] => Dict(:lifted_var_ref::Expr, :key::Expr, :lifted_constr_ref::Expr)
		[:(x[1]),:(x[2])] => [:refpod] = :(x[5]) <==> x[5] = x[1] * x[2]
		lifted_constr_ref is the lift constraint Expr
"""
function expr_t2y(terms::Array=[], t2y=Dict(), xdim::Int=0; kwargs...)

	if xdim > 0
		t2y[:xdim] = xdim
	end
	yidx = t2y[:xdim] + length(t2y)
	for i in terms
		if length(i) > 1 && !(i in keys(t2y)) && !(reverse(i) in keys(t2y)) #Perform lifting on multi-linear terms
			liftvarref = Expr(:ref, :x, yidx)
            liftconstr = Expr(:call, :(==), liftvarref, Expr(:call, :*, i[1], i[2])) # Bilinear
			t2y[i] = Dict(:lifted_var_ref=>liftvarref, :ref=>i, :lifted_constr_ref=>liftconstr)
			yidx += 1
		end
	end

	return t2y
end

"""
	Replace expression variable with lifted variable in-place
	Input: expression::Expr | t2y::Dict
	Output: recursive in-place operation
"""
function expr_lift(expr, t2y::Dict; kwargs...)

	@assert expr_mark(deepcopy(expr)) <= 2 # Not focusing on multi-linear terms
	expr_flatten(expr)	# Preproces the expression
	if expr.args[1] in [:+,:-]  || !isempty([i for i in 1:length(expr.args) if (isa(expr.args[i], Expr) && expr.args[i].head == :call)])
		for i in 2:length(expr.args) # Reserved 1 with operator
			if isa(expr.args[i], Expr) && expr.args[i].head == :call
				rep = expr_lift(expr.args[i], t2y) # Continue with sub-tree
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
		if haskey(t2y, refs)
			deleteat!(expr.args, idxs)
			push!(expr.args, t2y[refs][:lifted_var_ref]) # In-place Lift
		elseif haskey(t2y, reverse(refs)) # No duplicates
			deleteat!(expr.args, idxs)
			push!(expr.args, t2y[reverse(refs)][:lifted_var_ref]) # Place Lift
		end

	elseif expr.args[1] == :^  # with assumption :: square
		refs = [v for v in expr.args if (isa(v, Expr) && v.head == :ref)]
		power = [v for v in expr.args if (isa(v, Float64) || isa(v, Int64))]
		@assert length(refs) == 1  # Error handling
		@assert length(power) == 1
		@assert power[1] == 2
		refs = [refs;refs]
		@assert haskey(t2y, refs)
		expr = t2y[refs][:lifted_var_ref] # Lift
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
