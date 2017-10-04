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
