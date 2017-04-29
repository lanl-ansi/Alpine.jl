# This serves as a testing scripts for algorithmic functions
using MathProgBase, JuMP, Ipopt

include("../examples/bi1.jl")
include("../src/nlexpr.jl")
# include("src/formulator.jl")

m=pod_example_bi1()
d=JuMP.NLPEvaluator(m)
MathProgBase.initialize(d, [:ExprGraph])
ex=MathProgBase.constr_expr(d,2)
