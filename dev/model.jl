"""
	This is an experimental function that rebuilds the model with lifted variables.
	Later move this function to POD-based environment.
"""
function _pod_formulate_liftmodel(d, verbose=false)

	# !!!!!! Highly Experimental !!!!!! #

	# Collect all expressions
	expr_obj = deepcopy(MathProgBase.obj_expr(d))
	expr_cons = []
	for i in 1:length(d.constraints)
		push!(expr_cons, deepcopy(MathProgBase.constr_expr(d,i)))
	end

	TERMS, T2Y = pod_expr_rebuild(d, expr_obj)
	for i in 1:length(d.constraints)
		TERMS, T2Y = pod_expr_rebuild(d, expr_cons[i].args[2], TERMS, T2Y)
	end

    if verbose
        println("----------------------------------")
        for ml in keys(T2Y)
        	println("ML term $ml mapped to $(T2Y[ml][:podref])")
        end
        println("----------------------------------")
        println("OBJECTIVE => $(expr_obj)")
        println("----------------------------------")
        for i in 1:length(d.constraints)
        	println("CONS[$i] => $(expr_cons[i])")
        end
        println("----------------------------------")
    end

	pod = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0))) # Leave this as an option of the configuration
	podxdim = d.m.numCols + length(T2Y)
	@variable(pod, x[1:podxdim+1])

	for i in 1:d.m.numCols
		setcategory(x[i], d.m.colCat[i])
		setlowerbound(x[i], d.m.colLower[i])
		setupperbound(x[i], d.m.colUpper[i])
	end

	for i in keys(T2Y)
		JuMP.addNLconstraint(pod, T2Y[i][:link])
	end

	for i in expr_cons
		JuMP.addNLconstraint(pod, i)
	end

	JuMP.addNLconstraint(pod, Expr(:call, :(==), expr_obj, Expr(:ref, :x, (podxdim+1))))

	@objective(pod, d.m.objSense, x[end])
	if verbose
		print(pod)
	end

	pod.ext[:map] = T2Y

	return pod
end

"""
 	This is an experimental function that fetch variable informaiton.
	Later move this function POD-based environment.
"""
function pod_get_x_info(m, T2X)

	for i in keys(T2X)
		T2X[i][:lb] = getlowerbound(eval(T2X[i][:ref]))
		T2X[i][:ub] = getupperbound(eval(T2X[i][:ref]))
		T2X[i][:cate] = getcategory(eval[T2X[i][:ref]])
	end

	return T2X
end
