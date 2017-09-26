function ebd_sos1(L::Int, coding=ebd_binary)

	map = Dict()
	bCnt = Int(ceil(log(2,L)))
	cCnt = Int(2 * bCnt)
	for c in 1:cCnt
		map[c] = []
	end

	lseq = collect(0:(L-1))
	for l in 1:(L)	#Loop through all partitions in binary sense
		binCode = coding(lseq[l], bCnt) # Convert to binary string
		for b in 1:bCnt
			binCode[b] == '0' ? push!(map[b], lseq[l]+1) : push!(map[bCnt + b], lseq[l]+1)
		end
	end

	return map
end

# Make sure two things :: correct setup & compatible mapping function
function ebd_sos2(λCnt::Int, coding=ebd_gray)

	# Initialization
	map = Dict()
	logCnt = Int(ceil(log(2, λCnt-1)))
	for i in 1:logCnt*2 map[i]=[] end
	code_seq = [coding(i, logCnt) for i in 0:max(1,(2^logCnt-1))]

	# @show λCnt, logCnt, (logCnt^2-1), code_seq

	# Compatible Checking
	!is_compatible_encoding(code_seq) && error("Encodign method is not SOS-2 compatible...")

	# Mapping SOS-2 :: λ + λ <= x  &&  λ + λ <=  1-x
	# |x| = L : Flip  1 -> bCnt : bCnt + 1 -> bCnt + bCnt
	for l in 1:logCnt
		for j in 0:(λCnt-1)
			prod([(l in ebd_σ(code_seq[i])) for i in ebd_I(j, λCnt)]) && push!(map[l], j+1)
			prod([!(l in ebd_σ(code_seq[i])) for i in ebd_I(j, λCnt)]) && push!(map[l+logCnt], j+1)
		end
	end

	# @show map
	return map
end

function ebd_σ(b::String)
	sv = ebd_support_vec(b)
	return [i for i in 1:length(sv) if sv[i]]
end

function ebd_I(j, λCnt)
	return [i for i in collect(1:Int((λCnt-1))) if j in [i-1, i]] #Page 52
end

function is_compatible_encoding(code_seq::Vector)
	for i in 1:(length(code_seq)-1)
		sum(abs.(ebd_support_vec(code_seq[i])-ebd_support_vec(code_seq[i+1]))) != 1 && return false
	end
	return true
end

function ebd_support_vec(s::String)
   v = Vector{Bool}(length(s))
   for i in 1:length(s)
       s[i] == '1' ? v[i]=true : v[i]=false
   end
   return v
end

function ebd_binary(n, bCnt)
    return bin(n, bCnt)
end

# Given encoding bin
function ebd_gray(n, bCnt)
	code_decimal = xor(n, n >> 1)
	return bin(code_decimal, bCnt)
end
