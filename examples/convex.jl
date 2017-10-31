function convex_test(;verbose=false, solver=nothing, recognize=true)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
								   mip_solver=GurobiSolver(OutputFlag=1),
                                   discretization_ratio=4,
								   recognize_convex=recognize,
								   presolve_bound_tightening=false,
								   log_level=200))
	else
		m = Model(solver=solver)
	end

    @variable(m, 0<=x[1:5]<=2)

    @constraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                             # true
    @constraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                         # true
    @constraint(m, 3(x[1]x[1]) + 4*x[2]*x[2] <= -5)                             # false
    @constraint(m, 3(x[1]x[1]) + 4*x[2]^2 <= 10)                                # true
    @constraint(m, 3x[1]^2 + 4x[2]^2 + 6x[3]^2 <= 10)                           # true

	@NLconstraint(m, 3x[1]^0.5 + 4x[2]^0.5 + 5x[5]^0.5 <= 100)					# true | type-C
	@NLconstraint(m, -3x[1]^0.5 -4x[2]^0.5 >= -100)								# true | type-C

    @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] <= 25)                           # true

    @NLconstraint(m, (3*x[1]*x[1] + 4*x[2]*x[2]) <= 25)                         # true
    @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] - 25 <= 0)                       # true
    @NLconstraint(m, -3*x[1]*x[1] -4*x[2]*x[2] >= -25)                          # true
    @NLconstraint(m, 3*x[1]*x[1] + 5x[2]*x[2] <= 25)                            # true

    @NLconstraint(m, 4*x[1]^2 + 5x[2]^2 <= 25)                                  # Pass
    @NLconstraint(m, 3*x[1]*x[1] - 25 + 4*x[2]*x[2] <= 0)                       # false (unsupported when with @NLconstraint)
    @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[1] <= 25)                           # false
    @NLconstraint(m, 3*x[1]*x[1] + 16*x[2]^2 <= 40)                             # true
    @NLconstraint(m, 3*x[1]^2 + 16*x[2]^2 + 17 <= 16)                           # false

    @NLconstraint(m, 3*x[1]^3 + 16*x[2]^2 <= 20 - 20)                           # false
    @NLconstraint(m, 3*x[1]*x[1] + 4*x[2]*x[2] + 5*x[3]*x[3] + 6x[4]x[4] <= 15) # true
    @NLconstraint(m, 3x[1]x[1] + 4x[2]x[2] + 5x[3]^2 <= -15)                    # false
    @NLconstraint(m, 3x[1]^2 + 4x[2]^2 >= 15)                                   # false
    @NLconstraint(m, - 3x[1]^2 - 4x[3]^2 >= -15)                                # true

    @NLconstraint(m, 3x[1]^4 + 4x[2]^4 <= 200)                                  # true
    @NLconstraint(m, 3x[1]^4 + 4x[2]x[2]x[2]x[2] - 200 <= 0)                    # true
    @NLconstraint(m, 3x[1]^4 + 4x[2]^2*x[2]*x[2] <= 200)                        # true
    @NLconstraint(m, 3x[1]^4 + 4x[2]^3 <= 200)                                  # false
    @NLconstraint(m, 3x[1]^8 + 16*25*x[2]^8 - 30x[3]^8 <= 50)                   # false


    @objective(m, Max, x[1]^2+x[3]^2)

	if verbose
		print(m)
	end

	return m
end

function convex(;verbose=false, solver=nothing, recognize=true)

	return
end
