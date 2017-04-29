using MathProgBase, JuMP, Ipopt

include("../examples/bi2.jl")
include("../src/expr.jl")
include("../dev/model.jl")

bi_m=example_model(true)

bi_d = JuMP.NLPEvaluator(bi_m)
MathProgBase.initialize(bi_d, [:ExprGraph])

pod_m = _pod_formulate_liftmodel(bi_d, true)

solve(pod_m)
solve(bi_m)
