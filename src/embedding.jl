"""
	This function creates a mapping dictionary to link the right λ to the right bineary variables
	based on how many partitions are required and a given encoding method.
"""
function embedding_map(λCnt::Int, encoding::Any=ebd_gray, ibs::Bool=false)

	map = Dict()

	encoding = Alp.resolve_encoding_key(encoding)
	L = Int(ceil(log(2, λCnt-1)))
	for i in 1:L*2 map[i]=Set() end
	H = [encoding(i, L) for i in 0:max(1,(2^L-1))]
	map[:H_orig] = H
	map[:H] = [Alp.ebd_support_binary_vec(H[i]) for i in 1:length(H)] 		# Store the map
	map[:L] = L

	!is_compatible_encoding(H) && error("Encodign method is not SOS-2 compatible...")

	# Mapping SOS-2 :: Aλ <= x  &&  A'λ <=  1-x
	# |x| = L : Flip  1 -> bCnt : bCnt + 1 -> bCnt + bCnt
	J = Set(collect(1:(2^L+1)))
	J_s = Set(collect(1:λCnt))
	for l in 1:L
		if ibs
			L_constr = Set()
			R_constr = Set()
			for j in 0:2^L
				prod([(l in ebd_σ(H[i])) for i in ebd_I(j, 2^L+1)]) && push!(L_constr, j+1)
				prod([!(l in ebd_σ(H[i])) for i in ebd_I(j, 2^L+1)]) && push!(R_constr, j+1)
			end
			L_branch = intersect(J_s, setdiff(J, L_constr))
			R_branch = intersect(J_s, setdiff(J, R_constr))
			map[l] = setdiff(J_s, L_branch)
			map[l+L] = setdiff(J_s, R_branch)
		else
			for j in 0:(λCnt-1)
				prod([(l in ebd_σ(H[i])) for i in ebd_I(j, λCnt)]) && push!(map[l], j+1)
				prod([!(l in ebd_σ(H[i])) for i in ebd_I(j, λCnt)]) && push!(map[l+L], j+1)
			end
		end
	end

	return map
end

"""
	This function parse the encoding key string and return the built-in encoding function.
"""
function resolve_encoding_key(encoding::Any)
	isa(encoding, Function) && return encoding
	encoding == "default" && return ebd_gray
	error("Must specify a encoding method when using convhull_ebd formulation")
end

"""
	This function is the same σ() function described in Vielma and Nemhauser 2011.
"""
function ebd_σ(b::String)
	sv = Alp.ebd_support_bool_vec(b)
	return [i for i in 1:length(sv) if sv[i]]
end

"""
	This function is the same I() function described in Vielma and Nemhauser 2011.
"""
function ebd_I(j, λCnt)
	return [i for i in collect(1:Int((λCnt-1))) if j in [i-1, i]] #Page 52
end

"""
	This function is the same S() function described in Vielma and Nemhauser 2011.
"""
function ebd_S(i)
	(i <= 0) && error("Embedding utility S_[i] doesn't take i<=0")
	return [i, i-1]
end

"""
	This function checks whether the encoding method is compatible or not to obtain a valid mapping.
"""
function is_compatible_encoding(code_seq::Vector)
	for i in 1:(length(code_seq)-1)
		sum(abs.(Alp.ebd_support_bool_vec(code_seq[i])-Alp.ebd_support_bool_vec(code_seq[i+1]))) != 1 && return false
	end
	return true
end

"""
	This is a utility function that convert the binary string into bool vector
"""
function ebd_support_bool_vec(s::String)
   v = Vector{Bool}(undef, length(s))
   for i in 1:length(s)
       s[i] == '1' ? v[i]=true : v[i]=false
   end
   return v
end

"""
	This is a utility function that convert the binary string into 0/1 vector
"""
function ebd_support_binary_vec(s::String)
   v = Vector{Int}(undef, length(s))
   for i in 1:length(s)
       s[i] == '1' ? v[i]=1 : v[i]=0
   end
   return v
end

"""
	This is the function that translate the bounding constraints (α¹b⁰+α²b¹ <= x <= α¹b¹+α²b²)
	with log # of binary variables, i.e., generate these constraints using log # of binary variables.
"""
function ebd_link_xα(m::Optimizer, α::Vector, λCnt::Int, disc_vec::Vector, code_seq::Vector, var_idx::Int)

	lifters = Dict()
	exprs = Dict()
	L = Int(ceil(log(2, λCnt-1)))
	P = length(disc_vec) - 1

	# Expression expansion
	for i in 1:P
		code_vec = Alp.ebd_support_bool_vec(code_seq[i])
		lifters, exprs = Alp.ebd_link_expression(code_vec, lifters, exprs, i)
	end

	# Construct Variable Vector
	α_A = JuMP.@variable(m.model_mip, [1:length(keys(lifters))], lower_bound=0.0, upper_bound=1.0, base_name="αA$(var_idx)")
	for i in keys(lifters) # Build first-level evaluation
		Alp.binprod_relax(m.model_mip, α_A[lifters[i]-L], [α[j] for j in i])
	end

	α_R = [α; α_A] # Initialize/re-arrgange the variable sequence

	for i in 1:P # Start populating sub-expressions
		exprs[i][:expr] = JuMP.@expression(m.model_mip, sum(exprs[i][:coefs][j]*α_R[exprs[i][:vars][j]] for j in 1:exprs[i][:length] if exprs[i][:vars][j] != 0) + exprs[i][:vals])
	end

	# Contructing final constraints
	JuMP.@constraint(m.model_mip, _index_to_variable_ref(m.model_mip, var_idx) >= sum(exprs[j][:expr]*disc_vec[j] for j in 1:P))
	JuMP.@constraint(m.model_mip, _index_to_variable_ref(m.model_mip, var_idx) <= sum(exprs[j-1][:expr]*disc_vec[j] for j in 2:(P+1)))

	return
end

"""
	This is a utility function used in ebd_link_xα() that constructing the mapping that links partitions with
	combinations of binary variables.
"""
function ebd_link_expression(code::Vector, lift_dict::Dict, link_dict::Dict, p_idx::Int)

    L =  length(code)

    # Strickly mapping p_idx -> x
    code[1] ? coefs = Any[1] : coefs = Any[1, -1]
    code[1] ? vars = Any[1] : vars = Any[0, 1]
    prod(.!code) ? vals = 1.0 : vals = 0.0  # Marks the eventual expression is going have 1 or not

    for i in 2:L
        if code[i]
             for j in 1:length(coefs)
                (vars[j] == 0) ? vars[j] = i : vars[j] = [vars[j]; i] # Update vars
                @assert length(vars) == length(coefs)
            end
        else
            for j in 1:length(coefs)
                (vars[j] == 0) ? push!(vars, i) : push!(vars, [vars[j]; i]) # Update vars with variable multiplier (ADDING)
                (vars[j] == 0) ? push!(coefs, -1) : push!(coefs, coefs[j] * -1) # Update the coeffs with -1 (ADDING)
            end
        end
    end
    @assert length(coefs) == length(vars)

    # Maintaining the dictionary of lifted variables and swp
    for i in 1:length(vars)
        if isa(vars[i], Vector)
            !haskey(lift_dict, vars[i]) && (lift_dict[vars[i]] = L + length(keys(lift_dict)) + 1)
            vars[i] = lift_dict[vars[i]]
        end
    end

    # Maintain the expression dictionary
    link_dict[p_idx] = Dict(:p_idx=>p_idx, :coefs=>coefs, :vars=>vars, :vals=>vals, :length=>length(coefs))

    return lift_dict, link_dict
end

"""
	Built-in Encoding methods: binary encoding
	This is a compatible encoding
"""
ebd_binary(n, idx) = string(n, base=2, pad=idx) 

"""
	Built-in Encoding methods: gray encoding
	This is a compatible encoding
"""
function ebd_gray(n, idx)
	code_decimal = xor(n, n >> 1)
	return string(code_decimal, base=2, pad=idx)
end
