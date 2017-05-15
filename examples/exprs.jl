# Contains a basic model with various expressions for testing
using POD, JuMP, Ipopt, CPLEX, MathProgBase 

function example_exprs(verbose=false)

	m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0), mip_solver=CplexSolver()))
	@variable(m, px[i=1:6]>=1) # At some point if an initial value is given, keep them

	@NLconstraint(m, sum(3*px[i]^2 for i=1:4) >= 111)
	@NLconstraint(m,  - px[1]*px[2] + 4*5*px[3]*px[4] >= 222)
	@NLconstraint(m, -px[1]*px[2] <= 115)
	@NLconstraint(m, -px[1]*-px[2] >= 115)
	@NLconstraint(m, px[1]*-px[2] <= 115)
	@constraint(m, -px[1] + (-5) - 4<= 100)
	@NLconstraint(m, px[1]+ px[2]*px[3] >= 555) # => px[1] + x23 >= 555 && x23 == px[2]*px[3]
	@NLconstraint(m, px[1]^2 - 7*px[2]^2 + px[3]^2 + px[4] <= 6666)
	@NLconstraint(m, 13*px[1] - px[2] + 5*px[3]*6 + px[4] >= 77)

	@NLobjective(m, Min, 7*px[1]*6*px[4]*2+5+17+px[1]+px[2]+px[3]+8+3*5*px[1]^2*4)
	
	if verbose
		print(m)
	end

	return m
end
