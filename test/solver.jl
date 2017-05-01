using POD
using Ipopt
using JuMP
using CPLEX
initval = [1,5,5,1]
m = Model(solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0), log_level=1))
@variable(m, 1 <= x[i=1:4] <= 5, start=initval[i])
@constraint(m, sum(x) >= 0)
@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 25)
@NLconstraint(m, sum(x[i]^2 for i=1:4) == 40)
@NLobjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])
solve(m)

