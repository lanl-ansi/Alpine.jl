using JuMP, MathProgBase, Cbc, Ipopt, POD

function operator_b(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(OutputFlag=0),
								   presolve_bound_tightening=true,
								   presolve_bound_tightening_algo=2,
								   presolve_bt_output_tol=1e-1,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	@variable(m, x[1:4]>=0)
	@variable(m, y[1:3]<=0)
																# PARSING TARGETS #
	@constraint(m, x[1] + x[2] + -y[1] >= 1)					# None
	@constraint(m, 5x[1] + 8x[2] + 9y[2] <= 2)					# None
	@constraint(m, 4*x[1] + 5*x[2] + -3*y[1] >= 3)				# None
	@constraint(m, 3 <= x[1] + x[2] <= 4)						# None

	@NLconstraint(m, x[1]*x[2] + 3*x[2]*x[3] >= 5)				# x[1]*x[2] and 3*x[2]*x[3]
	@NLconstraint(m, x[1]*(x[2]+4) + 10*44*x[2]*5*x[3] <= 6)	# x[2]*x[3]
	@NLconstraint(m, x[1]*4 + x[2]*(x[2]+x[3]) >= 7)			# None
	@NLconstraint(m, 3*x[1]*x[2] + x[2]*-x[3] <= 8)				# x[1]*x[2]

	@NLconstraint(m, -1*x[1]^2 - x[2]^2 >= 9)					# x[1]^2 and x[2]^2
	@NLconstraint(m, x[1]^(1+1) + 5*x[2]^2 + 3*x[3]^3 <= 10)	# x[1]^2 and x[2]^2
	@NLconstraint(m, (x[1]+5)^2 >= 11)							# None
	@NLconstraint(m, x[1]^3 + x[3]^99 >= 12)					# None

	@NLconstraint(m, x[1]*x[2]*x[3] <= 13)						# x[1]*x[2]*x[3]
	@NLconstraint(m, (x[1]*x[2])*x[3] >= 14)					# x[1]*x[2]  ******
	@NLconstraint(m, x[1]*(x[2]*x[3]) <= 15)					# x[2]*x[3]  ******
	@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 16)					# x[1]*x[2]*x[3]*x[4]
	@NLconstraint(m, (x[1]*x[2])*(x[3]*x[4]) <= 17)				# x[1]*x[2] and x[3]*x[4]
	@NLconstraint(m, x[1]*(x[2]*x[3])*x[4] >= 18)				# x[2]*x[3]  ******
	@NLconstraint(m, (x[1]*x[2])*x[3]*x[4] <= 19)				# x[1]*x[2]  ******
	@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 20)				# x[1]*x[2]	 ******
	@NLconstraint(m, (x[1]*(x[2]*x[3])*x[4]) <= 21)				# x[2]*x[3]	 ******

	@NLconstraint(m, 4*5*6*x[1]*x[2]*x[3]*x[4] >= 22)			# x[1]*x[2]*x[3]*x[4]
	@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4] <= 23)			# None
	@NLconstraint(m, x[1]*x[1]*x[2]*x[2]*x[3]*x[3] >= 24)		# x[1]*x[1]*x[2]*x[2]*x[3]*x[3]
	@NLconstraint(m, (x[1]+1)*(x[2]+2)*(x[3]+3)*(x[4]+4) <= 25)	# None

	@NLconstraint(m, 50sin(x[1]) - 32*cos(x[2]) >= 26)			# sin(x[1]), cos(x[2])
	@NLconstraint(m, sin(x[1]+x[2]) - cos(x[2]-x[3]) + sin(-x[2]+x[3]) <= 27)
	@NLconstraint(m, sin(4*x[1]*x[2]) + cos(x[2]*-1*x[3]) >= 28)# x[1]*x[2] and x[2]*-1*x[3]  ******
	@NLconstraint(m, sin(x[1]*x[2]*x[3]) <= 29)					# x[1]*x[2]*x[3]

	if verbose
		print(m)
	end

	return m

end
