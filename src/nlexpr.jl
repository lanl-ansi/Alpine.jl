const operator_coefficient_base = Dict{Symbol, Float64}(
	:+ => 0.0,
	:- => 0.0,
	:/ => 1.0,
	:* => 1.0,
	:^ => 1.0
)

function pod_expr_rebuild(d, h, term::Array=[], t2y::Dict=Dict(); kwargs...)
	term, non = pod_expr_traversal(h, term)
	t2y = pod_expr_t2y(d.m.numCols, term, t2y)
	pod_expr_lift(h, t2y)
	return term, t2y
end

function pod_expr_t2y(xdim, T, t2y=Dict(); kwargs...)

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


function pod_expr_lift(h, t2y; kwargs...)

	@assert pod_expr_mark(deepcopy(h)) <= 2 # Not focusing on multi-linear termss
	pod_expr_flatten(h)
	if h.args[1] == :+  || !isempty([i for i in 1:length(h.args) if (isa(h.args[i], Expr) && h.args[i].head == :call)])
		for i in 2:length(h.args)
			if isa(h.args[i], Expr) && h.args[i].head == :call
				rep = pod_expr_lift(h.args[i], t2y) # Extend the replacement
				h.args[i] = rep
			elseif isa(h.args[i], Float64) || h.args[i].head == :ref
				h.args[i] = h.args[i]
			end
		end

	elseif h.args[1] == :* # with some assumption
		idxs = [i for i in 1:length(h.args) if (isa(h.args[i], Expr) && h.args[i].head == :ref)]
		refs = h.args[idxs]
		if length(refs) == length(h.args) - 1
			push!(h.args, 1.0)
		end
		if haskey(t2y, refs)
			deleteat!(h.args, idxs)
			push!(h.args, t2y[refs][:podref])
		elseif haskey(t2y, reverse(refs))
			deleteat!(h.args, idxs)
			push!(h.args, t2y[reverse(refs)][:podref])
		end

	elseif h.args[1] == :^
		refs = [v for v in h.args if (isa(v, Expr) && v.head == :ref)]
		@assert length(refs) == 1
		refs = [refs;refs]
		@assert haskey(t2y, refs)
		h = t2y[refs][:podref]

	else
		error("Unspported operator for variable lifting $(h.args[1])")
	end

	return h
end

function pod_expr_flatten(h; kwargs...)

	flat = pod_expr_arrangeargs(h.args)
	if isa(flat, Float64)
		return flat
	else
		h.args = flat
	end

	for i in 2:length(h.args)
		if isa(h.args[i], Expr) && h.args[i].head == :call
			h.args[i] = pod_expr_flatten(h.args[i])
		end
	end

	h.args = pod_expr_arrangeargs(h.args)

	return h
end

function pod_expr_arrangeargs(args; kwargs...)

	reverter = Dict(true=>:-, false=>:+)

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
		val = operator_coefficient_base[args[1]]
		exprs = []
		refs = []
		valinit = false
		refinit = false
		callinit = false
		for i in 2:length(args)
			if isa(args[i],Float64)
				valinit += (i==2)
				val = eval(args[1])(val, eval(reverter[(i==2) && (args[1]==:-)])(args[i]))
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
		if val == operator_coefficient_base[args[1]]
			val = []
		end
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

function pod_expr_mark(h; kwargs...)

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
				p = pod_expr_dim_enquiry(h)
			elseif h.args[1] in [:/]
				error("Does not support :/ and :^ yet.")
			end
		elseif isa(h.args[i], Expr) && h.args[i].head == :call
			if h.args[1] in [:+,:-]
				innerp = pod_expr_mark(h.args[i])
				p = max(p, innerp)
				h.args[i].typ = innerp
			elseif h.args[1] in [:*]
				innerp = pod_expr_mark(h.args[i])
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

function pod_expr_dim_enquiry(h)
	@assert h.args[1] == :^
	@assert length(h.args) == 3
	if isa(h.args[2], Float64)
		return h.args[2]
	else
		return h.args[3]
	end
end

function pod_expr_islinear(h)
	if pod_expr_mark(h) > 1
		return False
	else
		return True
	end
end

function pod_expr_traversal(h, T=[], buffer=[]; kwargs...)

	if isa(h, Float64) || isa(h, Int64)
		return T, buffer
	elseif h == :+ || h == :-
		if !(buffer in T) && !isempty(buffer)
			push!(T,buffer)
		end
		buffer = []
		return T, buffer
	elseif h == :*
		return T, buffer
	elseif h == :^
		return T, buffer
	elseif h.head == :ref
		push!(buffer, h)
		return T, buffer
	end

	for t in 1:length(h.args)
		T, buffer = pod_expr_traversal(h.args[t], T, buffer)
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
