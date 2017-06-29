using JuMP, MathProgBase
using Ipopt, Cbc
using JLD
using POD

m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
						   mip_solver=CbcSolver(),
						   presolve_perform_bound_tightening=true,
						   presolve_bound_tightening_algo=2,
						   presolve_bt_output_tolerance=1e-1,
						   log_level=1))

@variable(m, x[1:4]>=0)

@NLconstraint(m, x[1]^2 >= 1)
@NLconstraint(m, x[1]*x[2] <= 1)
@NLconstraint(m, x[1]^2 + x[2]*x[3] <= 1)

@NLconstraint(m, x[1] * (x[2]*x[3]) >= 1)
@NLconstraint(m, x[1]^2 * (x[2]^2 * x[3]^2) <= 1)

@NLconstraint(m, (x[1] * x[2]) * x[3] >= 1)
@NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

@NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
@NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2 <= 1)
@NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
@NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 1)
@NLconstraint(m, ((x[1]^2*x[2])*x[3])*x[4] <= 1)
@NLconstraint(m, ((x[1]*x[2]^2)*x[3])*x[4] >= 1)
@NLconstraint(m, ((x[1]*x[2])*x[3]^2)*x[4] <= 1)
@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4]^2 >= 1)
@NLconstraint(m, ((x[1]^2*x[2]^2)*x[3]^2)*x[4]^2 <=1)

@NLconstraint(m, x[1]*(x[2]*(x[3]*x[4])) >= 1)
@NLconstraint(m, x[1]^2*(x[2]*(x[3]*x[4])) <= 1)
@NLconstraint(m, x[1]*(x[2]^2*(x[3]*x[4])) >= 1)
@NLconstraint(m, x[1]*(x[2]*(x[3]^2*x[4])) <= 1)
@NLconstraint(m, x[1]*(x[2]*(x[3]*x[4]^2)) >= 1)
@NLconstraint(m, x[1]^2*(x[2]^2*(x[3]^2*x[4]^2)) <= 1)

@NLconstraint(m, x[1]*x[2]*x[3] >= 1)
@NLconstraint(m, x[1]^2*x[2]*x[3] >= 1)
@NLconstraint(m, x[1]*x[2]^2*x[3] >= 1)
@NLconstraint(m, x[1]*x[2]*x[3]^2 >= 1)
@NLconstraint(m, x[1]^2*x[2]^2*x[3] >= 1)
@NLconstraint(m, x[1]*x[2]^2*x[3]^2 >= 1)
@NLconstraint(m, x[1]^2*x[2]*x[3]^2 >= 1)
@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2 >= 1)

@NLconstraint(m, (x[1]*x[2])*(x[3]*x[4]) >= 1)
@NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]) >= 1)
@NLconstraint(m, (x[1]*x[2]^2)*(x[3]*x[4]) >= 1)
@NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]) >= 1)
@NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]^2) >= 1)
@NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]^2) >= 1)
@NLconstraint(m, (x[1]^2*x[2])*(x[3]^2*x[4]) >= 1)

@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 1)
@NLconstraint(m, x[1]^2*x[2]*x[3]*x[4] >= 1)
@NLconstraint(m, x[1]*x[2]^2*x[3]*x[4] >= 1)
@NLconstraint(m, x[1]*x[2]*x[3]^2*x[4]^2 >= 1)
@NLconstraint(m, x[1]*x[2]^2*x[3]^2*x[4] >= 1)
@NLconstraint(m, x[1]^2*x[2]*x[3]*x[4]^2 >= 1)
@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4]^2 >= 1)

@NLconstraint(m, (x[1]*x[2]*x[3])*x[4] >= 1)
@NLconstraint(m, (x[1]^2*x[2]*x[3])*x[4] >= 1)
@NLconstraint(m, (x[1]*x[2]^2*x[3])*x[4] >= 1)
@NLconstraint(m, (x[1]*x[2]*x[3]^2)*x[4] >= 1)
@NLconstraint(m, (x[1]*x[2]*x[3])*x[4]^2 >= 1)
@NLconstraint(m, (x[1]^2*x[2]^2*x[3]^2)*x[4]^2 >= 1)

@NLconstraint(m, x[1]*(x[2]*x[3])*x[4] >= 1)
@NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4] >= 1)
@NLconstraint(m, x[1]*(x[2]^2*x[3])*x[4] >= 1)
@NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4] >= 1)
@NLconstraint(m, x[1]*(x[2]*x[3])*x[4]^2 >= 1)
@NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4]^2 >= 1)
@NLconstraint(m, x[1]*(x[2]^2*x[3]^2)*x[4] >= 1)
@NLconstraint(m, x[1]^2*(x[2]^2*x[3])*x[4] >= 1)
@NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4]^2 >= 1)

@NLconstraint(m, x[1]*(x[2]*x[3]*x[4]) >= 1)
@NLconstraint(m, x[1]^2*(x[2]*x[3]*x[4]) >= 1)
@NLconstraint(m, x[1]*(x[2]^2*x[3]*x[4]) >= 1)
@NLconstraint(m, x[1]*(x[2]*x[3]^2*x[4]) >= 1)
@NLconstraint(m, x[1]*(x[2]*x[3]*x[4]^2) >= 1)
@NLconstraint(m, x[1]^2*(x[2]*x[3]^2*x[4]) >= 1)
@NLconstraint(m, x[1]^2*(x[2]^2*x[3]^2*x[4]^2) >= 1)

@NLconstraint(m, x[1]*x[2]*(x[3]*x[4]) >= 1)
@NLconstraint(m, x[1]^2*x[2]*(x[3]*x[4]) >= 1)
@NLconstraint(m, x[1]*x[2]^2*(x[3]*x[4]) >= 1)
@NLconstraint(m, x[1]*x[2]*(x[3]^2*x[4]) >= 1)
@NLconstraint(m, x[1]*x[2]*(x[3]*x[4]^2) >= 1)
@NLconstraint(m, x[1]^2*x[2]^2*(x[3]^2*x[4]^2) >= 1)

@NLconstraint(m, (x[1]*x[2])*x[3]*x[4] >= 1)
@NLconstraint(m, (x[1]^2*x[2])*x[3]*x[4] >= 1)
@NLconstraint(m, (x[1]*x[2]^2)*x[3]*x[4] >= 1)
@NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4] >= 1)
@NLconstraint(m, (x[1]*x[2])*x[3]*x[4]^2 >= 1)
@NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4]^2 >= 1)
@NLconstraint(m, (x[1]^2*x[2]^2)*x[3]^2*x[4]^2 >= 1)

JuMP.build(m)
