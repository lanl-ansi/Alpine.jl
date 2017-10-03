function embedding_sos1(L::Int, encoding::Any=ebd_binary)

	map = Dict()
	bCnt = Int(ceil(log(2,L)))
	cCnt = Int(2 * bCnt)
	for c in 1:cCnt
		map[c] = []
	end

	encoding = resolve_encoding_key(encoding)

	lseq = collect(0:(L-1))
	for l in 1:(L)	#Loop through all partitions in binary sense
		binCode = encoding(lseq[l], bCnt) # Convert to binary string
		for b in 1:bCnt
			binCode[b] == '0' ? push!(map[b], lseq[l]+1) : push!(map[bCnt + b], lseq[l]+1)
		end
	end

	return map
end

function embedding_sos1b(L::Int, encoding::Any=ebd_binary)

	map = Dict()
	bCnt = Int(ceil(log(2,L)))
	cCnt = Int(2 * bCnt)
	for c in 1:cCnt
		map[c] = []
	end

	encoding = resolve_encoding_key(encoding)

	lseq = collect(0:(L-1))
	entry_cnt = 0
	for l in 1:(L)	#Loop through all partitions in binary sense
		entry_cnt += 1
		λ_idx = [l,l+1]
		binCode = encoding(lseq[l], bCnt) # Convert to binary string
		J = ebd_σ(binCode)
		map[entry_cnt] = Dict(:L=>λ_idx, :Y=>J)
	end

	return map
end

# Make sure two things :: correct setup & compatible mapping function
function embedding_sos2(λCnt::Int, encoding::Any=ebd_gray)

	# Initialization
	encoding = resolve_encoding_key(encoding)
	map = Dict()
	logCnt = Int(ceil(log(2, λCnt-1)))
	for i in 1:logCnt*2 map[i]=[] end
	code_seq = [encoding(i, logCnt) for i in 0:max(1,(2^logCnt-1))]

	# @show λCnt, logCnt, (logCnt^2-1), code_seq
	# Compatible Checking
	!is_compatible_encoding(code_seq) && error("Encodign method is not SOS-2 compatible...")

	# Mapping SOS-2 :: Aλ <= x  &&  A'λ <=  1-x
	# |x| = L : Flip  1 -> bCnt : bCnt + 1 -> bCnt + bCnt
	# Theorem 3
	for l in 1:logCnt
		for j in 0:(λCnt-1)
			prod([(l in ebd_σ(code_seq[i])) for i in ebd_I(j, λCnt)]) && push!(map[l], j+1)
			prod([!(l in ebd_σ(code_seq[i])) for i in ebd_I(j, λCnt)]) && push!(map[l+logCnt], j+1)
		end
	end

	return map
end

function resolve_encoding_key(encoding::Any)

	isa(encoding, Function) && return encoding

	encoding == "default" && return ebd_gray
	encoding == "gray" && return ebd_gray
	encoding == "binary" && return ebd_binary
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
