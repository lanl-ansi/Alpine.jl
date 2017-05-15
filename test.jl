include("examples/exprs.jl")
m=example_exprs()
solve(m)
print(m)
print(m.internalModel.basic_model_mip)

for k in 1:m.internalModel.num_constr_orig
    ex = m.internalModel.lifted_constr_expr_mip[k]
    println(ex)
    println(POD.expr_to_affine(ex))
end