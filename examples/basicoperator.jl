using JuMP
using MathProgBase
using JLD
using POD

m = Model()

@variable(m, x[1:4]>=0)

@NLconstraint(m, x[1]^2 >= 1)  					# Basic monomial x[5]=x[1]^2

@NLconstraint(m, x[1]*x[2] <= 1)				#
@NLconstraint(m, x[1]^2 * x[2]^2 <= 1)

@NLconstraint(m, x[1]*(x[2]*x[3]) >= 1)
@NLconstraint(m, x[1]^2*(x[2]^2 * x[3]^3) <= 1)

@NLconstraint(m, (x[1]*x[2]) * x[3] >= 1)
@NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

@NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
@NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2) <= 1)
@NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
@NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

@NLconstraint(m, (((x[1]*x[2])*x[3])*x[4] >= 1)
@NLconstraint(m, (((x[1]^2*x[2])*x[3])*x[4] <= 1)
@NLconstraint(m, (((x[1]*x[2]^2)*x[3])*x[4] >= 1)
@NLconstraint(m, (((x[1]*x[2])*x[3]^2)*x[4] <= 1)
@NLconstraint(m, (((x[1]*x[2])*x[3])*x[4]^2 >= 1)
@NLconstraint(m, (((x[1]^2*x[2]^2)*x[3]^2)*x[4]^2) <=1)

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


d = JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:ExprGraph])

exprs = []
for i in 1:1
	ex = MathProgBase.constr_expr(d, i)
	@show i,ex
	push!(exprs, ex)
end
