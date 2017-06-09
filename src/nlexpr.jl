"""
	Populate the dictionary that contains the information of all the non-linear terms
"""
function populate_nonlinear_info(m::PODNonlinearModel; kwargs...)

	m.lifted_obj_expr_mip = deepcopy(m.obj_expr_orig)
 	for i in 1:m.num_constr_orig
		push!(m.lifted_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
	end

	m.nonlinear_info[:xdim] = m.num_var_orig
	terms = []
	expr_resolve_sign(m.lifted_obj_expr_mip)
	analyze_expr(m.lifted_obj_expr_mip, terms, m.nonlinear_info)
	for i in 1:m.num_constr_orig
		expr_resolve_sign(m.lifted_constr_expr_mip[i])
		analyze_expr(m.lifted_constr_expr_mip[i].args[2], terms, m.nonlinear_info)
	end

	delete!(m.nonlinear_info,:xdim)
	return m
end

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
		buffer = []  # Clear buffer
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

	if expr.args[1] in [:/]
		error("Unsupported expression operator $expr.args[1] during expression traversal")
	end

	return terms, buffer
end

"""
	This is wrapper generating an affine data structure from a lifted linear expression
	The lifted expressions are stored in the POD model as lifted_obj_expr_mip and lifted_constr_expr_mip
"""
function populate_lifted_affine(m::PODNonlinearModel; kwargs...)

	# Populate the affine data structure for the objective
	m.lifted_obj_aff_mip = expr_to_affine(m.lifted_obj_expr_mip)

	# Populate the affine data structure for the constraints
	for i in 1:m.num_constr_orig
		push!(m.lifted_constr_aff_mip, expr_to_affine(m.lifted_constr_expr_mip[i]))
	end

end

"""
	This function takes a constraint/objective expression and converts it into a affine expression data structure
	Use the function to traverse linear expressions traverse_expr_to_affine()
"""
function expr_to_affine(expr)

	# The input should follow :(<=, LHS, RHS)
	affdict = Dict()
	if expr.args[1] in [:(==), :(>=), :(<=)] # For a constraint expression
		@assert isa(expr.args[3], Float64) || isa(expr.args[3], Int)
		@assert isa(expr.args[2], Expr)
		# non are buffer spaces, not used anywhere
		lhscoeff, lhsvars, rhs, non, non = traverse_expr_to_affine(expr.args[2])
		rhs = -rhs + expr.args[3]
		affdict[:sense] = expr.args[1]
	else # For an objective expression
		lhscoeff, lhsvars, rhs, non, non = traverse_expr_to_affine(expr)
		affdict[:sense] = nothing
	end

	affdict[:coefs] = lhscoeff
	affdict[:vars]	= lhsvars
	affdict[:rhs]	= rhs
	@assert length(affdict[:coefs]) == length(affdict[:vars])
	affdict[:cnt]	= length(affdict[:coefs])

	return affdict
end

"""
	This function traverse a left hand side tree to collect affine terms
"""
function traverse_expr_to_affine(expr, lhscoeffs=[], lhsvars=[], rhs=0.0, bufferVal=0.0, bufferVar=nothing, level=0; kwargs...)

	reversor = Dict(true => -1.0, false => 1.0)

	if isa(expr, Float64) || isa(expr, Int) # Capture any coefficients or right hand side
		bufferVal = expr
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:+, :-]
		if bufferVal != 0.0 && bufferVar != nothing
			push!(lhscoeffs, bufferVal)
			push!(lhsvars, bufferVar)
			bufferVal = 0.0
			bufferVar = nothing
		end
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:*, :(<=), :(==), :(>=)]
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:/, :^]
		error("Unsupported operators $expr")
	elseif expr.head == :ref
		bufferVar = expr
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	end

	for i in 1:length(expr.args)
		lhscoeff, lhsvars, rhs, bufferVal, bufferVar=
			traverse_expr_to_affine(expr.args[i], lhscoeffs, lhsvars, rhs, bufferVal, bufferVar, level+1)
		if expr.args[1] in [:+, :-]  # Term segmentation [:-, :+], see this and close the previous term
			if bufferVal != 0.0 && bufferVar != nothing
				push!(lhscoeffs, reversor[(expr.args[1]==:- && i<=2)]*eval(expr.args[1])(bufferVal))
				push!(lhsvars, bufferVar)
				bufferVal = 0.0
				bufferVar = nothing
			end
			if bufferVal != 0.0 && bufferVar == nothing
				rhs = eval(expr.args[1])(rhs, bufferVal)
				bufferVal = 0.0
			end
			if bufferVal == 0.0 && bufferVar != nothing && expr.args[1] == :+
				push!(lhscoeffs, 1.0)
				push!(lhsvars, bufferVar)
				bufferVar = nothing
			end
			if bufferVal == 0.0 && bufferVar != nothing && expr.args[1] == :-
				push!(lhscoeffs, reversor[(i<=2)]*(-1.0))
				push!(lhsvars, bufferVar)
				bufferVar = nothing
			end
		elseif expr.args[1] in [:(<=), :(==), :(>=)]
			rhs = expr.args[end]
		end
	end

	if level == 0
		if bufferVal != 0.0 && bufferVar != nothing
			push!(lhscoeffs, bufferVal)
			push!(lhsvars, bufferVar)
			bufferVal = 0.0
			bufferVar = nothing
		end
	end

	return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
end

"""
	Formulate the pod mip model lifted constraints
"""
function populate_lifted_expr(m::PODNonlinearModel; kwargs...)

	expr_lift(m.lifted_obj_expr_mip, m.nonlinear_info)
	for i in 1:m.num_constr_orig
		expr_lift(m.lifted_constr_expr_mip[i].args[2], m.nonlinear_info)
	end

	return m
end

"""
	Dereferencing function: splice and dereference (not used anywhere)
"""
function expr_dereferencing(expr, m)

	for i in 2:length(expr.args)
		if isa(expr.args[i], Float64)
			k = 0
		elseif expr.args[i].head == :ref
			@assert isa(expr.args[i].args[2], Int)
			expr.args[i] = Variable(m, expr.args[i].args[2])
		elseif expr.args[i].head == :call
			expr_dereferencing(expr.args[i], m)
		else
			error("expr_dereferencing :: Unexpected term in expression tree.")
		end
	end

end

"""
	Analyze a single expression and create lifted terms and in-place lifted expression

	Input: expression | existing terms already parsed | nonlinear_info
	Output: Updated terms array and nonlinear_info (not returned, updated data structures)
	Algorithm:
		1) traverse the expression for term detection using traverse_expr
		2) lift the original expression in place

"""
function analyze_expr(expr, terms::Array=[], nonlinear_info::Dict=Dict(); kwargs...)

	terms, buffer = traverse_expr(expr, terms)
	add_info_using_terms(terms, nonlinear_info)

	return
end


"""
	Adding the lifted variable map using the non-linear terms

	Input: array of non-linear terms | nonlinear_info
	Output: updated nonlinear_info (no explicit return)

	Mapping : term => Dict(:lifted_var_ref::Expr, :ref::Expr, :lifted_constr_ref::Expr, :monomial_status::Bool)

"""
function add_info_using_terms(terms::Array=[], nonlinear_info=Dict(); kwargs...)

	@assert haskey(nonlinear_info, :xdim)
	yidx = nonlinear_info[:xdim] + length(nonlinear_info)
	# Logic only works for bi-linear terms
	for term in terms
		if length(term) > 1 && !(term in keys(nonlinear_info)) && !(reverse(term) in keys(nonlinear_info)) # Check for duplicates
			lifted_var_ref = Expr(:ref, :x, yidx)
            lifted_constr_ref = Expr(:call, :(==), lifted_var_ref, Expr(:call, :*, term[1], term[2]))
            monomial_status = false
            if term[1] == term[2]
                monomial_status = true
            end
            nonlinear_info[term] = Dict(:lifted_var_ref => lifted_var_ref,
                                        :ref => term,
                                        :lifted_constr_ref => lifted_constr_ref,
                                        :monomial_status => monomial_status,
										:cos_status => false,
										:sin_status => false)
			yidx += 1
		end
	end

	return
end

# from here on this code was written by sitew and has not been verified.

"""
	Replace expression variable with lifted variable in-place
	Input: expression::Expr | nonlinear_info::Dict
	Output: recursive in-place operation
"""
function expr_lift(expr, nonlinear_info::Dict; kwargs...)

	@assert expr_mark(deepcopy(expr)) <= 2 # Not focusing on multi-linear terms
	expr_flatten(expr)	# Preproces the expression
	if expr.args[1] in [:+,:-]  || !isempty([i for i in 1:length(expr.args) if (isa(expr.args[i], Expr) && expr.args[i].head == :call)])
		for i in 2:length(expr.args) # Reserved 1 with operator
			if isa(expr.args[i], Expr) && expr.args[i].head == :call
				rep = expr_lift(expr.args[i], nonlinear_info) # Continue with sub-tree
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
		if haskey(nonlinear_info, refs)
			deleteat!(expr.args, idxs)
			push!(expr.args, nonlinear_info[refs][:lifted_var_ref]) # In-place Lift
		elseif haskey(nonlinear_info, reverse(refs)) # No duplicates
			deleteat!(expr.args, idxs)
			push!(expr.args, nonlinear_info[reverse(refs)][:lifted_var_ref]) # Place Lift
		end

	elseif expr.args[1] == :^  # with assumption :: square
		refs = [v for v in expr.args if (isa(v, Expr) && v.head == :ref)]
		power = [v for v in expr.args if (isa(v, Float64) || isa(v, Int64))]
		@assert length(refs) == 1  # Error handling
		@assert length(power) == 1
		@assert power[1] == 2
		refs = [refs;refs]
		@assert haskey(nonlinear_info, refs)
		expr = nonlinear_info[refs][:lifted_var_ref] # Lift
	else
		error("Unspported operator for variable lifting $(expr.args[1])")
	end

	return expr
end


"""
	This function is trates the sign in expressions trees by cleaning the following structure:
		args
		::+ || ::-
		::Call || ::ref
	By separating the structure with some dummy treatments
"""
function expr_resolve_sign(expr, level=0; kwargs...)

	resolver = Dict(:- => -1, :+ => 1)
	for i in 2:length(expr.args)
		if !isa(expr.args[i], Float64) && !isa(expr.args[i], Int) 								# Skip the coefficients
			if length(expr.args[i].args) == 2 && expr.args[i].head == :call
				if expr.args[i].args[1] == :- 													# Treatment for :(-x), replace with :(x*-1)
					if expr.args[1] in [:*] 													# Only perform the treatment when connected with *
						push!(expr.args, resolver[expr.args[i].args[1]])
						expr.args[i] = expr.args[i].args[2]
					elseif expr.args[1] in [:+, :-, :(==), :(<=), :(>=)]
						expr.args[i] = expr.args[i]
					else
						error("Unexpected operator $(expr.args[1]) during resolving sign")
					end
				elseif expr.args[i].args[1] == :+ 												# Treatement for :(+x) replace with :x
					expr.args[i] = expr.args[i].args[2]
				end
			elseif expr.args[i].head == :call
				expr_resolve_sign(expr.args[i], level+1)
			end
		end
	end

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

	if level > 0  #Root level process, no additional operators
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
			if isa(args[i],Float64) || isa(args[i], Int)
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
		elseif isa(expr.args[i],Float64) || isa(expr.args[i],Int)
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
