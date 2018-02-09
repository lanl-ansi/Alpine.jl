"""
    process_expr(expr; kwargs...)

High-level warpper for processing expression with sub-tree operators
"""
function process_expr(m::PODNonlinearModel)

    expr_initialization(m)      # S0 : initialize the space for parsing and analyzing
    expr_preprocess(m)          # S1 : pre-process the negative sign in expressions
    expr_parsing(m)             # S2 : parsing the expressions for nonlinear information
    expr_conversion(m)          # S3 : convert lifted(linear) expressions into affine function
    expr_finalized(m)           # S4 : finalize process by extracting some measurements

    return
end

"""
	STEP 1: initialize the expression/ space
"""
function expr_initialization(m::PODNonlinearModel)

	# 0 : deepcopy data into mip lifted expr place holders
	m.bounding_obj_expr_mip = deepcopy(m.obj_expr_orig)
	m.bounding_obj_mip = Dict()

	for i in 1:m.num_constr_orig
		push!(m.bounding_constr_expr_mip, deepcopy(m.constr_expr_orig[i]))
		push!(m.bounding_constr_mip, Dict())
	end

	return
end

"""
	STEP 2: preprocess expression for trivial sub-trees and nasty pieces for easier later process
"""
function expr_preprocess(m::PODNonlinearModel)

	expr_resolve_const(m.bounding_obj_expr_mip)
	expr_resolve_sign(m.bounding_obj_expr_mip)
	expr_flatten(m.bounding_obj_expr_mip)
	for i in 1:m.num_constr_orig
		expr_resolve_const(m.bounding_constr_expr_mip[i])
		expr_resolve_sign(m.bounding_constr_expr_mip[i])
		expr_flatten(m.bounding_constr_expr_mip[i].args[2])
	end

	return
end

"""
	STEP 3: parse expression for patterns on either the generic level or term level
"""
function expr_parsing(m::PODNonlinearModel)

	is_strucural = expr_constr_parsing(m.bounding_obj_expr_mip, m)
	if !is_strucural
		m.bounding_obj_expr_mip = expr_term_parsing(m.bounding_obj_expr_mip, 0, m)
		m.obj_structure = :generic_linear
	end
	(m.loglevel > 99) && println("[OBJ] $(m.obj_expr_orig)")

	for i in 1:m.num_constr_orig
		is_strucural = expr_constr_parsing(m.bounding_constr_expr_mip[i], m, i)
		if !is_strucural
			m.bounding_constr_expr_mip[i] = expr_term_parsing(m.bounding_constr_expr_mip[i], i, m)
			m.constr_structure[i] = :generic_linear
		end
		(m.loglevel > 99) && println("[CONSTR] $(m.constr_expr_orig[i])")
	end

	return
end

"""
	STEP 4: convert the parsed expressions into affine-based function that can be used for adding JuMP constraints
"""
function expr_conversion(m::PODNonlinearModel)

	if m.obj_structure == :generic_linear
		m.bounding_obj_mip = expr_linear_to_affine(m.bounding_obj_expr_mip)
		m.obj_structure = :affine
	end
	m.loglevel > 99 && println("type :: ", m.obj_structure)
	m.loglevel > 99 && println("lifted ::", m.bounding_obj_expr_mip)
	m.loglevel > 99 && println("coeffs ::", m.bounding_obj_mip[:coefs])
	m.loglevel > 99 && println("vars ::", m.bounding_obj_mip[:vars])
	m.loglevel > 99 && println("sense ::", m.bounding_obj_mip[:sense])
	m.loglevel > 99 && println("rhs ::", m.bounding_obj_mip[:rhs])
	m.loglevel > 99 && println("----------------")


	for i in 1:m.num_constr_orig
		if m.constr_structure[i] == :generic_linear
			m.bounding_constr_mip[i] = expr_linear_to_affine(m.bounding_constr_expr_mip[i])
			m.constr_structure[i] = :affine
		end
		m.loglevel > 99 && println("type :: ", m.constr_structure[i])
		m.loglevel > 99 && println("lifted ::", m.bounding_constr_expr_mip[i])
		m.loglevel > 99 && println("coeffs ::", m.bounding_constr_mip[i][:coefs])
		m.loglevel > 99 && println("vars ::", m.bounding_constr_mip[i][:vars])
		m.loglevel > 99 && println("sense ::", m.bounding_constr_mip[i][:sense])
		m.loglevel > 99 && println("rhs ::", m.bounding_constr_mip[i][:rhs])
		m.loglevel > 99 && println("----------------")
	end

	return
end


"""
	STEP 5: collect measurements and information as needed for handly operations in the algorithm section
"""
function expr_finalized(m::PODNonlinearModel)

	collect_nonlinear_vars(m)
	m.candidate_disc_vars = sort(m.candidate_disc_vars)
    m.num_var_linear_mip = length(m.linear_terms)
	m.num_var_nonlinear_mip = length(m.nonlinear_terms)
	m.num_constr_convex = length([i for i in m.constr_structure if i == :convex])

	return m
end

function collect_nonlinear_vars(m::PODNonlinearModel)

    # Walk through all nonlinear terms
	for i in keys(m.nonlinear_terms)
        m.nonlinear_terms[i][:discvar_collector](m, i)
    end

	return
end

function expr_strip_const(expr, subs=[], rhs=0.0)

    exhaust_const = [!(expr.args[1] in [:+, :-]) || !(isa(expr.args[i], Float64) || isa(expr.args[i], Int)) for i in 2:length(expr.args)]
    if prod(exhaust_const)
        push!(subs, expr)
        return subs, rhs
    end

    for i in 2:length(expr.args)
        if (isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol))
            (expr.args[1] == :+) && (rhs -= expr.args[i])
            (expr.args[1] == :-) && ((i == 2) ? rhs -= expr.args[i] : rhs += expr.args[i])
        elseif expr.args[i].head == :ref
            continue
        elseif expr.args[i].head == :call
            subs, rhs = expr_strip_const(expr.args[i], subs, rhs)
        end
    end
    return subs, rhs
end

"""
    This utility function builds a constraint reference with by repating one
    operator with a vector variable references.
    For example,
    input => y, x[1,2,3], :+
    output => (y = x[1] + x[2] + x[3])::Expr
"""
function build_constr_block(y_idx::Int, var_idxs::Vector, operator::Symbol)

    expr_l_block = Expr(:call, :(==), parse("x[$(y_idx)]"))
    expr_r_block = Expr(:call, operator)
    for j in 1:length(var_idxs)
        push!(expr_r_block.args, parse("x[$(var_idxs[j])]"))
    end
    push!(expr_l_block.args, expr_r_block)
    return expr_l_block
end

"""
    expr_constr_parsing(expr, m::PODNonlinearModel)

Recognize structural constraints.
"""
function expr_constr_parsing(expr, m::PODNonlinearModel, idx::Int=0)

    # First process user-defined structures in-cases of over-ride
    for i in 1:length(m.constr_patterns)
        is_strucural = eval(m.constr_patterns[i])(expr, m, idx)
        return
    end

    # Recognize built-in special structural pattern
    if m.recognize_convex
        is_convex = resolve_convex_constr(expr, m, idx)
        is_convex && return true
    end

    # More patterns goes here

    return false
end

function expr_is_axn(expr, scalar=1.0, var_idxs=[], power=[]; N=nothing)

    # println("inner recursive input ", expr)
    !(expr.args[1] in [:*,:^]) && return nothing, nothing, nothing         # Limited Area

    if expr.args[1] == :*
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                scalar *= expr.args[i]
            elseif (expr.args[i].head == :ref)
                @assert isa(expr.args[i].args[2], Int)
                push!(var_idxs, expr.args[i])
                push!(power, 1)
            elseif (expr.args[i].head == :call)
                scalar, var_idxs, power = expr_is_axn(expr.args[i], scalar, var_idxs, power)
            end
        end
    elseif expr.args[1] == :^
        for i in 2:length(expr.args)
            if isa(expr.args[i], Float64) || isa(expr.args[i], Int)
                push!(power, expr.args[i])
                continue
            elseif (expr.args[i].head == :ref)
                push!(var_idxs, expr.args[i])
                continue
            elseif (expr.args[i].head == :call)
                return nothing, nothing, nothing
            end
        end
    end

    # If the user wants a specific N
    !(N == nothing) && !(length(var_idxs) == N) && return nothing, nothing, nothing
    (var_idxs == nothing) && (scalar == nothing) && (power == nothing) && return nothing, nothing, nothing # Unrecognized sub-structure

    # println("inner recursive ", scalar, " ", var_idxs)
    @assert length(var_idxs) == length(power)
    return scalar, var_idxs, power
end

"""
	This function takes a constraint/objective expression and converts it into a affine expression data structure
	Use the function to traverse linear expressions traverse_expr_linear_to_affine()
"""
function expr_linear_to_affine(expr)

	# The input should follow :(<=, LHS, RHS)
	affdict = Dict()
	if expr.args[1] in [:(==), :(>=), :(<=)] # For a constraint expression
		@assert isa(expr.args[3], Float64) || isa(expr.args[3], Int)
		@assert isa(expr.args[2], Expr)
		# non are buffer spaces, not used anywhere
		lhscoeff, lhsvars, rhs, non, non = traverse_expr_linear_to_affine(expr.args[2])
		rhs = -rhs + expr.args[3]
		affdict[:sense] = expr.args[1]
	elseif expr.head == :ref  # For single variable objective expression
		lhscoeff = [1.0]
		lhsvars = [expr]
		rhs = 0
		affdict[:sense] = nothing
	else # For an objective expression
		lhscoeff, lhsvars, rhs, non, non = traverse_expr_linear_to_affine(expr)
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

	traverse_expr_linear_to_affine(expr, lhscoeffs=[], lhsvars=[], rhs=0.0, bufferVal=0.0, bufferVar=nothing, sign=1.0, level=0)

This function traverse a left hand side tree to collect affine terms.
Updated status : possible to handle (x-(x+y(t-z))) cases where signs are handled properly
"""
function traverse_expr_linear_to_affine(expr, lhscoeffs=[], lhsvars=[], rhs=0.0, bufferVal=0.0, bufferVar=nothing, sign=1.0, coef=1.0, level=0)

	# @show expr, coef, bufferVal

	reversor = Dict(true => -1.0, false => 1.0)
	function sign_convertor(subexpr, pos)
		if length(subexpr.args) == 2 && subexpr.args[1] == :-
			return -1.0
	 	elseif length(subexpr.args) > 2 && subexpr.args[1] == :- && pos > 2
			return -1.0
		end
		return 1.0
	end

	if isa(expr, Float64) || isa(expr, Int) # Capture any coefficients or right hand side
		(bufferVal > 0.0) ? bufferVal *= expr : bufferVal = expr * coef
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:+, :-]    # TODO: what is this condition?
		if bufferVal != 0.0 && bufferVar != nothing
			push!(lhscoeffs, bufferVal)
			push!(lhsvars, bufferVar)
			bufferVal = 0.0
			bufferVar = nothing
		end
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:*]
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:(<=), :(==), :(>=)]
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	elseif expr in [:/, :^]
		error("Unsupported operators $expr, it is suppose to be affine function")
	elseif expr.head == :ref
		bufferVar = expr
		return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
	end

	# HOTPATCH : Special Structure Recognization
	start_pos = 1
	if (expr.args[1] == :*) && (length(expr.args) == 3)
		if (isa(expr.args[2], Float64) || isa(expr.args[2], Int)) && (expr.args[3].head == :call)
			(coef != 0.0) ? coef = expr.args[2] * coef : coef = expr.args[2]	# Patch
			# coef = expr.args[2]
			start_pos = 3
			warn("Speicial expression structure detected [*, coef, :call, ...]. Currently handling using a beta fix...", once=true)
		end
	end

	for i in start_pos:length(expr.args)
		lhscoeff, lhsvars, rhs, bufferVal, bufferVar = traverse_expr_linear_to_affine(expr.args[i], lhscoeffs, lhsvars, rhs, bufferVal, bufferVar, sign*sign_convertor(expr, i), coef, level+1)
		if expr.args[1] in [:+, :-]  # Term segmentation [:-, :+], see this and wrap-up the current (linear) term
			if bufferVal != 0.0 && bufferVar != nothing  # (sign) * (coef) * (var) => linear term
				push!(lhscoeffs, sign*sign_convertor(expr, i)*bufferVal)
				push!(lhsvars, bufferVar)
				bufferVal = 0.0
				bufferVar = nothing
			end
			if bufferVal != 0.0 && bufferVar == nothing  # (sign) * (coef) => right-hand-side term
				rhs += sign*sign_convertor(expr, i)*bufferVal
				bufferVal = 0.0
			end
			if bufferVal == 0.0 && bufferVar != nothing && expr.args[1] == :+
				push!(lhscoeffs, sign*1.0*coef)
				push!(lhsvars, bufferVar)
				bufferVar = nothing
			end
			if bufferVal == 0.0 && bufferVar != nothing && expr.args[1] == :-
				push!(lhscoeffs, sign*sign_convertor(expr, i)*coef)
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
						error("Unexpected operator $(expr.args[i]) during resolving sign")
					end
				elseif expr.args[i].args[1] == :+ 												# Treatement for :(+x) replace with :x
					expr.args[i] = expr.args[i].args[2]
				end
			elseif expr.args[i].head == :call
				expr_resolve_sign(expr.args[i], level+1)
			end
		end
	end

	return
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

	if level > 0  #Root level process, no additional processes
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

	if args[1] in [:^, :sin, :cos]
		return args
	elseif args[1] in [:/]
		warn("Partially supported operator $(args[1]) in $(args). Trying to resolve the expression...")
		args = expr_resolve_divdenominator(args)
	elseif args[1] in [:*, :+, :-]
		# Such generic operators can be handled
	else
		error("Unkown/unspported operator $(args[1])")
	end

	if args[1] in [:*, :+, :-]
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
		if args[1] in [:*]
			return[args[1];val;refs;exprs]
		elseif args[1] in [:+, :-]
			if Bool(valinit)
				return [args[1];val;refs;exprs]
			elseif Bool(refinit)
				return [args[1];refs;val;exprs]
			elseif Bool(callinit)
				return [args[1];exprs;refs;val]
			else
				error("Unexpected condition encountered...")
			end
		end
	else
		error("Unsupported expression arguments $args")
	end

	return
end

"""
	Check if a sub-tree is pure constant or not
"""
function expr_resolve_const(expr)

	for i in 1:length(expr.args)
		if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
			continue
		elseif expr.args[i].head == :call
			(expr_isconst(expr.args[i])) && (expr.args[i] = eval(expr.args[i]))
			if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
				continue
			else
				expr_resolve_const(expr.args[i])
			end
		end
	end

	return
end

"""
	Check if a division's denominator is pure numerical and replace it with a multiply
"""
function expr_resolve_divdenominator(args)

	@assert args[1] == :/
	@assert length(args) == 3
	resolvable = true

	for i in 3:length(args)
		resolvable *= expr_isconst(args[i])
	end

	if resolvable
		for i in 3:length(args)
			args[i]=1/eval(args[i])
		end
		args[1] = :*
	end

	return args
end

"""
	Check if a sub-tree(:call) is totally composed of constant values
"""
function expr_isconst(expr)

	(isa(expr, Float64) || isa(expr, Int) || isa(expr, Symbol)) && return true

	const_tree = true
	for i in 1:length(expr.args)
		if isa(expr.args[i], Float64) || isa(expr.args[i], Int) || isa(expr.args[i], Symbol)
			continue
		elseif expr.args[i].head == :call
			const_tree *= expr_isconst(expr.args[i])
		elseif expr.args[i].head == :ref
			return false
		end
	end

	return const_tree
end
