include("examples/nlp1.jl")
using JLD
m = nlp1()
JuMP.build(m)
imp = load("examples/exprs.jld")
exprs = imp["exprs"]
for i in 1:length(exprs)
	@show i, exprs[i]
end

POD.expr_podding(exprs[5], m.internalModel)
POD.expr_podding(exprs[13], m.internalModel)
