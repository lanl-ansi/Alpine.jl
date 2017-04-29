# This is a function under development
function pod_formulate_liftmodel(d, verbose=false)

	# Collect all expressions
	expr_obj = MathProgBase.obj_expr(d)
	expr_cons = []
	for i in 1:length(d.constraints)
		push!(expr_cons, MathProgBase.constr_expr(d,i))
	end

	# Analyze Formulation
	TERMS, T2Y = pod_expr_rebuild(expr_obj)
	for i in 1:length(d.constraints)
		TERMS, T2Y = pod_expr_rebuild(expr_cons[i].args[2], TERMS, T2Y)
	end

	if verbose
		println("\n----------------------------------")
		for ml in keys(T2Y)
			println("ML term $ml mapped to $(T2Y[ml][:podref])")
		end
		println("\n----------------------------------")
		println("OBJECTIVE => $(expr_obj)")
		println("\n----------------------------------")
		for i in 1:length(d.constraints)
			println("CONS[$i] => $(expr_cons[i])")
		end
		println("\n----------------------------------")
	end

    #POD Model Constrauction
	pod = Model(solver=IpoptSolver()) # Leave this as an option of the configuration
	@variable(pod, x[1:(d.m.numCols+length(T2Y))])

	# Construct Lifting Links
	for i in keys(T2Y)
		println(typeof(i))
		@constraint(pod, T2Y[i][:link])
	end

	# Reconstruct Linear Model
	for i in expr_cons
		println(typeof(i))
		@constraint(pod, i)
	end

	@objective(pod, m.objSense, eval(expr_obj))
	if verbose
		print(pod)
	end
	return pod
end

function pod_get_x_info(m, T2X)
	for i in keys(T2X)
		T2X[i][:lb] = getlowerbound(eval(T2X[i][:ref]))
		T2X[i][:ub] = getupperbound(eval(T2X[i][:ref]))
		T2X[i][:cate] = getcategory(eval[T2X[i][:ref]])
	end
	return T2X
end
