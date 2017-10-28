# Make sure two things :: correct setup & compatible mapping function
function embedding_map(λCnt::Int, encoding::Any=ebd_gray, ibs::Bool=false)

	map = Dict()

	encoding = resolve_encoding_key(encoding)
	L = Int(ceil(log(2, λCnt-1)))
	for i in 1:L*2 map[i]=Set() end
	H = [encoding(i, L) for i in 0:max(1,(2^L-1))]
	map[:H_orig] = H
	map[:H] = [ebd_support_binary_vec(H[i]) for i in 1:length(H)] 		# Store the map
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

function resolve_encoding_key(encoding::Any)
	isa(encoding, Function) && return encoding
	encoding == "default" && return ebd_gray
	encoding == "gray" && return ebd_gray
	encoding == "binary" && return ebd_binary
	error("Must specify a encoding method when using embedding formulation")
end


function ebd_σ(b::String)
	sv = ebd_support_bool_vec(b)
	return [i for i in 1:length(sv) if sv[i]]
end

function ebd_I(j, λCnt)
	return [i for i in collect(1:Int((λCnt-1))) if j in [i-1, i]] #Page 52
end

function ebd_S(i)
	(i <= 0) && error("Embedding utility S_[i] doesn't take i<=0")
	return [i, i-1]
end

function is_compatible_encoding(code_seq::Vector)
	for i in 1:(length(code_seq)-1)
		sum(abs.(ebd_support_bool_vec(code_seq[i])-ebd_support_bool_vec(code_seq[i+1]))) != 1 && return false
	end
	return true
end

function ebd_support_bool_vec(s::String)
   v = Vector{Bool}(length(s))
   for i in 1:length(s)
       s[i] == '1' ? v[i]=true : v[i]=false
   end
   return v
end

function ebd_support_binary_vec(s::String)
   v = Vector{Int}(length(s))
   for i in 1:length(s)
       s[i] == '1' ? v[i]=1 : v[i]=0
   end
   return v
end

# Version 1 implementation #
function ebd_link_xα_v1(m::PODNonlinearModel, α::Vector, λCnt::Int, disc_vec::Vector, code_seq::Vector, var_idx::Int)

	var_refs = Dict()
	L = Int(ceil(log(2, λCnt-1)))
	P = length(disc_vec) - 1
	α_C = @variable(m.model_mip, [1:L], lowerbound=0.0, upperbound=1.0, basename="αC$(var_idx)")
	if L > 2
		α_S = Dict()
		for i in 1:P
			α_S[i] = @variable(m.model_mip, [2:L], lowerbound=0.0, upperbound=1.0, basename="αS$(var_idx)")
		end
	end
	α_L = @variable(m.model_mip, [1:P], lowerbound=0.0, upperbound=1.0, basename="αL$(var_idx)")

	# Linking (1 - x) = y
	@constraint(m.model_mip, [i in 1:L], α_C[i] == 1 - α[i])

    #            S in 2:L
	#            |
	#        ____|
	#        |   |
	# Form ((x * x) * x) ... * x) * x  --- i in 1:length(code_seq)
	#      |<-        L             |
	# Regulated x with α_C
	for i in 1:P
		code_vec = ebd_support_bool_vec(code_seq[i])
		if L == 2
			mccormick_bin(m.model_mip, α_L[i], code_vec[1] ? α[1] : α_C[1], code_vec[2] ? α[2] : α_C[2])
		else
			if !haskey(var_refs, (code_vec[1] ? α[1].col : α_C[1].col, code_vec[2] ? α[2].col : α_C[2].col))
				mccormick_bin(m.model_mip, α_S[i][2], code_vec[1] ? α[1] : α_C[1], code_vec[2] ? α[2] : α_C[2])
				var_refs[(code_vec[1] ? α[1].col : α_C[1].col, code_vec[2] ? α[2].col : α_C[2].col)] = α_S[i][2]
			else
				@show "duplicating A, key =$((code_vec[1] ? α[1].col : α_C[1].col, code_vec[2] ? α[2].col : α_C[2].col))"
				α_S[i][2] = var_refs[(code_vec[1] ? α[1].col : α_C[1].col, code_vec[2] ? α[2].col : α_C[2].col)]
			end
			for j in 3:L
				if !haskey(var_refs, (α_S[i][j-1].col, code_vec[j] ? α[j].col : α_C[j].col))
					mccormick_bin(m.model_mip, α_S[i][j], α_S[i][j-1], code_vec[j] ? α[j] : α_C[j])
					var_refs[(α_S[i][j-1].col, code_vec[j] ? α[j].col : α_C[j].col)] = α_S[i][j]
				else
					@show "duplicating A, key =$((α_S[i][j-1].col, code_vec[j] ? α[j].col : α_C[j].col))"
					α_S[i][j] = var_refs[(α_S[i][j-1].col, code_vec[j] ? α[j].col : α_C[j].col)]
				end
			end
			if !haskey(var_refs, (α_S[i][L-1].col, α_S[i][L].col))
				mccormick_bin(m.model_mip, α_L[i], α_S[i][L-1], α_S[i][L])
				var_refs[(α_S[i][L-1].col, α_S[i][L].col)] = α_L[i]
			else
				α_L[i] = var_refs[(α_S[i][L-1].col, α_S[i][L].col)]
			end
		end
	end
	@constraint(m.model_mip, Variable(m.model_mip, var_idx) >= sum(α_L[j]*disc_vec[j] for j in 1:P)) # Add x = f(α) for regulating the domains
	@constraint(m.model_mip, Variable(m.model_mip, var_idx) <= sum(α_L[j-1]*disc_vec[j] for j in 2:(P+1)))
	return
end


# Version 2
function ebd_link_xα_v2(m::PODNonlinearModel, α::Vector, λCnt::Int, disc_vec::Vector, code_seq::Vector, var_idx::Int)

	var_refs = Dict()
	L = Int(ceil(log(2, λCnt-1)))
	P = length(disc_vec) - 1
	α_C = @variable(m.model_mip, [1:L], lowerbound=0.0, upperbound=1.0, basename="αC$(var_idx)")
	α_L = @variable(m.model_mip, [1:P], lowerbound=0.0, upperbound=1.0, basename="αL$(var_idx)")

	# Construct the (1 - x) = y
	@constraint(m.model_mip, [i in 1:L], α_C[i] == 1 - α[i])

	# Regulated x with α_C
	for i in 1:P
		code_vec = ebd_support_bool_vec(code_seq[i])
		binprod_relax(m.model_mip, α_L[i], [code_vec[j] ? α[j] : α_C[j] for j in 1:L])
	end

	@constraint(m.model_mip, Variable(m.model_mip, var_idx) >= sum(α_L[j]*disc_vec[j] for j in 1:P)) # Add x = f(α) for regulating the domains
	@constraint(m.model_mip, Variable(m.model_mip, var_idx) <= sum(α_L[j-1]*disc_vec[j] for j in 2:(P+1)))

	return
end

# Version 3
function ebd_link_xα_v3(m::PODNonlinearModel, α::Vector, λCnt::Int, disc_vec::Vector, code_seq::Vector, var_idx::Int)

	lifters = Dict()
	exprs = Dict()
	L = Int(ceil(log(2, λCnt-1)))
	P = length(disc_vec) - 1

	# Expression expansion
	for i in 1:P
		code_vec = ebd_support_bool_vec(code_seq[i])
		lifters, exprs = ebd_link_expression(code_vec, lifters, exprs, i)
	end

	# Construct Variable Vector
	α_A = @variable(m.model_mip, [1:length(keys(lifters))], lowerbound=0.0, upperbound=1.0, basename="αA$(var_idx)")
	for i in keys(lifters) # Build first-level evaluation
		binprod_relax(m.model_mip, α_A[lifters[i]-L], [α[j] for j in i])
	end

	α_R = [α, α_A;] # Iintialize Rearrgange the variable sequence

	for i in 1:P # Start populating sub-expressions
		exprs[i][:expr] = @expression(m.model_mip, sum(exprs[i][:coefs][j]*α_R[exprs[i][:vars][j]] for j in 1:exprs[i][:length] if exprs[i][:vars][j] != 0) + exprs[i][:vals])
	end

	# Contructing final constraints
	@constraint(m.model_mip, Variable(m.model_mip, var_idx) >= sum(exprs[j][:expr]*disc_vec[j] for j in 1:P))
	@constraint(m.model_mip, Variable(m.model_mip, var_idx) <= sum(exprs[j-1][:expr]*disc_vec[j] for j in 2:(P+1)))

	return
end

function ebd_link_expression(code::Vector, lift_dict::Dict, link_dict::Dict, p_idx::Int)

    L =  length(code)

    # Strickly mapping p_idx -> x
    code[1] ? coefs = Any[1] : coefs = Any[1, -1]
    code[1] ? vars = Any[1] : vars = Any[0, 1]
    prod(.!code) ? vals = 1.0 : vals = 0.0  # Marks the eventual expression is going have 1 or not

    for i in 2:L
        if code[i]
             for j in 1:length(coefs)
                (vars[j] == 0) ? vars[j] = i : vars[j] = [vars[j], i;] # Update vars
                @assert length(vars) == length(coefs)
            end
        else
            for j in 1:length(coefs)
                (vars[j] == 0) ? push!(vars, i) : push!(vars, [vars[j], i;]) # Update vars with variable multiplier (ADDING)
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

# ============================= #
#   Built-in Encoding methods   #
# ============================= #
function ebd_binary(n, idx)
    return bin(n, idx)
end

function ebd_gray(n, idx)
	code_decimal = xor(n, n >> 1)
	return bin(code_decimal, idx)
end
