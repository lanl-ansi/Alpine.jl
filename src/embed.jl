function p2B_sos1(L::Int, coding=encode_binary)

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
function p2B_sos2(λCnt::Int, coding=encoding_gray)

	# Initialization
	logCnt = Int(ceil(log(2,λCnt-1)))
	map = Dict()
	code_seq = [coding(i) for i in 0:λCnt]

	# Compatible Checking
	!is_compatible_encoding(code_seq) && error("Encodign method is not SOS-2 compatible...")

	# Mapping SOS-2 :: λ + λ <= x  &&  λ + λ <=  1-x
	# |x| = L : Flip  1 -> bCnt : bCnt + 1 -> bCnt + bCnt
	for 1 in 1:logCnt
		# Positve selection
		# 0 selection
	end

	return map
end

# Output J vector for selection
function ebd_J(l)
	return
end

function ebd_B()
	return
end

function ebd_I(j, λCnt)
	return [i for i in collect(1:(λCnt-1)) if j in [i-1, i]] #Page 52
end

function is_compatible_encoding(code_seq::Vector)

	for i in 1:(length(code_seq)-1)
		sum(code_seq[i] .!= code_seq[i+1]) != 1 return false
	end

	return true
end

function encode_binary(n, bCnt)
    return bin(n, bCnt)
end

# Given encoding bin
function encode_gray(n, bCnt)
	code_decimal = xor(n, n >> 1)
	return bin(code_decimal, bCnt)
end


function p2B_convhull(L::Int)
    return
end
