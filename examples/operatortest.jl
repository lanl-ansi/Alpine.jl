include("examples/advanceoperator.jl")
JuMP.build(m)
@show exprs[5]
POD.expr_parsing(exprs[5], m.internalModel)
@show exprs[6]
POD.expr_parsing(exprs[6], m.internalModel)
@show exprs[9]
POD.expr_parsing(exprs[9], m.internalModel)
@show exprs[13]
POD.expr_parsing(exprs[13], m.internalModel)
@show exprs[17]
POD.expr_parsing(exprs[17], m.internalModel)
