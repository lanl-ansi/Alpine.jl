# Make sure two things :: correct setup & compatible mapping function
function embedding_map(λCnt::Int, encoding::Any=ebd_gray, ibs::Bool=false)

	map = Dict()

	encoding = resolve_encoding_key(encoding)
	L = Int(ceil(log(2, λCnt-1)))
	for i in 1:L*2 map[i]=Set() end
	H = [encoding(i, L) for i in 0:max(1,(2^L-1))]
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


function ebd_link_xα(m::PODNonlinearModel, α::Vector, λCnt::Int, disc_vec::Int, code_seq::Vector, var_idx::Int)

	L = Int(ceil(log(2, λCnt-1)))
	P = length(disc_vec) - 1
	α_C = @variable(m.model_mip, [1:L], lowerbound=0.0, upperbound=1.0, basename="αC$(var_idx)")
	L > 2 && (α_S = @variable(m.model_mip, [2:L], lowerbound=0.0, upperbound=1.0, basename="αS$(var_idx)"))
	α_L = @variable(m.model_mip, [1:(length(disc_vec)-1)], lowerbound=0.0, upperbound=1.0, basename="αL$(var_idxs)")

	# Linking (1 - x) = y
	@variable(m.model_mip, [for i in 1:L], α_C[i] == 1 - α[i])

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
			mccormick_bin(m.model_mip, α_L[i], code_vec[1] ? α[1] : α_S[1], code_vec[1] ? α[2] : α_S[2])
		else
			mccormick_bin(m.model_mip, α_S[2], code_vec[1] ? α[1] : α_S[1], code_vec[2] ? α[2] : α_S[2])
			for j in 3:L
				mccormick_bin(m.model_mip, α_S[j], α_S[j-1], code_vec[j] ? α[j] : α_S[j])
			end
			mccormick_bin(m.model_mip, α_L[i], α_S[L-1], α_S[L])
		end
	end

	@constraint(m.model_mip, Variable(m.model_mip, var_idx) >= sum(α_L[j]*disc_vec[j] for j in 1:P)) # Add x = f(α) for regulating the domains
	@constraint(m.model_mip, Variable(m.model_mip, var_idx) <= sum(α_L[j-1]*diss_vec[j] for j in 2:(P+1)))

	return
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
