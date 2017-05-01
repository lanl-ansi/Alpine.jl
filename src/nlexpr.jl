function pod_expr_traversal(h, T::Array=[], buffer::Array=[]; kwargs...)

	if isa(h, Float64) || isa(h, Int64)
		return T, buffer  # Continue to collect
	elseif h == :+ || h == :-
		if !(buffer in T) && !isempty(buffer)
			push!(T,buffer)
		end
		buffer = []
		return T, buffer  # Term collecing
	elseif h == :*
		return T, buffer  # Continue to collect
	elseif h == :^
		return T, buffer # Continue to collect
	elseif h == :/
		error("Unspported operator $(h)")
	elseif h.head == :ref
		push!(buffer, h)
		return T, buffer
	end

	for t in 1:length(h.args)
		T, buffer = pod_expr_traversal(h.args[t], T, buffer)
		if h.args[1] == :+ || h.args[1] == :-
			if !(buffer in T) && !isempty(buffer)
				push!(T, buffer)
			end
			buffer = []
		end
	end

	# Post Process on special operators
	if h.args[1] == :^
		push!(buffer, buffer[1])
		if !(buffer in T) && !isempty(buffer)
			push!(T, buffer)
		end
		buffer = []
	end

	return T, buffer
end
