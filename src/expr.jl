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
function pod_expr_rebuild(d::JuMP.NLPEvaluator, h, term::Array=[], t2y::Dict=Dict(); kwargs...)

	options = Dict(kwargs)
	term, non = _pod_expr_traversal(h, term)
	t2y = _pod_expr_t2y(d.m.numCols, term, t2y)
	_pod_expr_lift(h, t2y)

	return term, t2y
end

"""
	Lift variable mapping :: initialize or continue build variable mapping

	Input: numCols::Int | lifted terms::Array | t2y mapping dictionary::Dict
	Output: updated mapping dictionary::Dict

	Naming : x[1:numCols] > [x;lifted x] > lifted x still use name x
	Mapping : [:Expr, :Expr] => Dict(:podref::Expr, :key::Expr, :link::Expr)
		[:(x[1]),:(x[2])] => [:refpod] = :(x[5]) <==> x[5] = x[1] * x[2]
		link is the lift constraint Expr
"""
function _pod_expr_t2y(xdim::Int, T::Array=[], t2y=Dict(); kwargs...)

	@assert xdim > 0
	yidx = xdim + length(t2y) + 1
	for i in T
		if length(i) > 1 && !(i in keys(t2y)) && !(reverse(i) in keys(t2y)) #Perform lifting on multi-linear terms
			liftvarref = Expr(:ref, :x, yidx)
            liftconstr = Expr(:call, :(==), liftvarref, Expr(:call, :*, i[1], i[2])) # Bilinear
			t2y[i] = Dict(:podref=>liftvarref, :ref=>i, :link=>liftconstr)
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
function _pod_expr_lift(h, t2y::Dict; kwargs...)

	@assert _pod_expr_mark(deepcopy(h)) <= 2 # Not focusing on multi-linear terms
	_pod_expr_flatten(h)	# Preproces the expression

	if h.args[1] in [:+,:-]  || !isempty([i for i in 1:length(h.args) if (isa(h.args[i], Expr) && h.args[i].head == :call)])
		for i in 2:length(h.args) # Reserved 1 with operator
			if isa(h.args[i], Expr) && h.args[i].head == :call
				rep = _pod_expr_lift(h.args[i], t2y) # Continue with sub-tree
				h.args[i] = rep
			elseif isa(h.args[i], Float64) || h.args[i].head == :ref
				h.args[i] = h.args[i]
			end
		end

	elseif h.args[1] == :* # with assumption :: bilinear, limited () usage
		idxs = [i for i in 1:length(h.args) if (isa(h.args[i], Expr) && h.args[i].head == :ref)]
		refs = h.args[idxs]
		if length(refs) == length(h.args) - 1  # (*)->[(x),(x)] ==> (*)->[(1),(y)]
			push!(h.args, 1.0)
		end
		if haskey(t2y, refs)
			deleteat!(h.args, idxs)
			push!(h.args, t2y[refs][:podref]) # In-place Lift
		elseif haskey(t2y, reverse(refs)) # No duplicates
			deleteat!(h.args, idxs)
			push!(h.args, t2y[reverse(refs)][:podref]) # Place Lift
		end

	elseif h.args[1] == :^  # with assumption :: square
		refs = [v for v in h.args if (isa(v, Expr) && v.head == :ref)]
		power = [v for v in h.args if (isa(v, Float64) || isa(v, Int64))]
		@assert length(refs) == 1  # Error handling
		@assert length(power) == 1
		@assert power[1] == 2
		refs = [refs;refs]
		@assert haskey(t2y, refs)
		h = t2y[refs][:podref] # Lift
	else
		error("Unspported operator for variable lifting $(h.args[1])")
	end

	return h
end

"""
	Recursivly pre-process expression by treating the coefficients
"""
function _pod_expr_flatten(h; kwargs...)

	flat = _pod_expr_arrangeargs(h.args)
	if isa(flat, Float64)
		return flat
	else
		h.args = flat
	end

	for i in 2:length(h.args)
		if isa(h.args[i], Expr) && h.args[i].head == :call
			h.args[i] = _pod_expr_flatten(h.args[i])
		end
	end

	h.args = _pod_expr_arrangeargs(h.args)

	return h
end

"""
	Re-arrange childrens by their type.
	Only consider 3 types: coefficients(Int64, Float64), :ref, :calls
	Sequence arrangement is dependent on operator
"""
function _pod_expr_arrangeargs(args::Array; kwargs...)

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
function _pod_expr_mark(h; kwargs...)

	p = 0
	for i in 2:length(h.args)
		if isa(h.args[i], Expr) && h.args[i].head == :ref
			if h.args[1] in [:+, :-]
				p = max(p, 1)
				h.args[i].typ = 1
			elseif h.args[1] in [:*]
				p = +(p, 1)
				h.args[i].typ = 1
			elseif h.args[1] in [:^] # Consider leaf when encounter :^
				p = _pod_expr_dim_enquiry(h)
			elseif h.args[1] in [:/]
				error("Does not support :/ and :^ yet.")
			end
		elseif isa(h.args[i], Expr) && h.args[i].head == :call
			if h.args[1] in [:+,:-]
				innerp = _pod_expr_mark(h.args[i])
				p = max(p, innerp)
				h.args[i].typ = innerp
			elseif h.args[1] in [:*]
				innerp = _pod_expr_mark(h.args[i])
				p = +(p, innerp)
				h.args[i].typ = innerp
			elseif h.args[1] in [:/, :^]
				error("Does not support :/ and :^ yet.")
			end
		elseif isa(h.args[i],Float64)
			p = p
		else
			error("Unexpected expression node encouter $(h.args[i])")
		end
	end

	return p
end

"""
	A special check when encouter :^. Cleaner function
"""
function _pod_expr_dim_enquiry(h)

	@assert h.args[1] == :^
	@assert length(h.args) == 3 # Don't take cray :call with :^
	if isa(h.args[2], Float64)
		return h.args[2]
	else
		return h.args[3]
	end

end

"""
	Check if a sub-tree is linear or not
"""
function pod_expr_islinear(h)
	if _pod_expr_mark(h) > 1
		return False
	else
		return True
	end
end

"""
	Tree traversal for bi-linear term detection and marking
	Input: expression subtree | terms::Array | buffer::Array
	Output: Found terms and an empty buffer

	Algorithm: DFS tree traversal to collect bi-linear terms using a buffer.
		When collected buffer encouter a breaker :(+) :(-), finish recognizing a term
		and collect it into T.

	NOTE: This function may change a lot when working on a multi-linear term version
"""
function _pod_expr_traversal(h, T::Array=[], buffer::Array=[]; kwargs...)

	if isa(h, Float64) || isa(h, Int64)
		return T, buffer  # Continue to collect
	elseif h == :+ || h == :-
		if !(buffer in T) && !isempty(buffer)
			push!(T,buffer)
		end
		buffer = []
		return T, buffer  # Term collecing
	elseif h == :*
		return T, buffer  # Continue to collect
	elseif h == :^
		return T, buffer # Continue to collect
	elseif h == :/
		error("Unspported operator $(h)")
	elseif h.head == :ref
		push!(buffer, h)
		return T, buffer
	end

	for t in 1:length(h.args)
		T, buffer = _pod_expr_traversal(h.args[t], T, buffer)
		if h.args[1] == :+ || h.args[1] == :-
			if !(buffer in T) && !isempty(buffer)
				push!(T, buffer)
			end
			buffer = []
		end
	end

	# Post Process on special operators
	if h.args[1] == :^
		push!(buffer, buffer[1])
		if !(buffer in T) && !isempty(buffer)
			push!(T, buffer)
		end
		buffer = []
	end

	return T, buffer
end
