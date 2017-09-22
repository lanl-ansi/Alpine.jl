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
