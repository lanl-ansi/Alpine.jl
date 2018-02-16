using Base.Test
using JuMP, MathProgBase
using Ipopt, Cbc, Pajarito
using GLPKMathProgInterface
using POD

poddir = Pkg.dir("POD")

examples = readdir(joinpath(Pkg.dir("POD"), "test", "examples"))
for i in examples
    include(joinpath(Pkg.dir("POD"), "test", "examples", i))
end

# Performe Tests
include("$(poddir)/test/solver.jl")
include("$(poddir)/test/expression.jl")
include("$(poddir)/test/algorithm.jl")
# include("$(poddir)/test/integer.jl")
include("$(poddir)/test/utility.jl")

# Test Generation Code commented
# for k in keys(m.internalModel.nonconvex_terms)
#     println("@test m.internalModel.nonconvex_terms[$(k)][:y_idx] == $(m.internalModel.nonconvex_terms[k][:y_idx])")
#     println("@test m.internalModel.nonconvex_terms[$(k)][:id] == $(m.internalModel.nonconvex_terms[k][:id])")
#     println("@test m.internalModel.nonconvex_terms[$(k)][:lifted_constr_ref] == :($(m.internalModel.nonconvex_terms[k][:lifted_constr_ref]))")
#     println("@test m.internalModel.nonconvex_terms[$(k)][:nonlinear_type] == :$(m.internalModel.nonconvex_terms[k][:nonlinear_type])")
#     println("@test m.internalModel.nonconvex_terms[$(k)][:y_type] == :$(m.internalModel.nonconvex_terms[k][:y_type])")
#     println("@test m.internalModel.nonconvex_terms[$(k)][:constr_id] == $(m.internalModel.nonconvex_terms[k][:constr_id])")
# end
#
#
# for (i,k) in enumerate(keys(m.internalModel.linear_terms))
#     println("lk$(i) = $(k)")
#     println("@test m.internalModel.linear_terms[lk$(i)][:y_idx] == $(m.internalModel.linear_terms[k][:y_idx])")
#     println("@test m.internalModel.linear_terms[lk$(i)][:id] == $(m.internalModel.linear_terms[k][:id])")
#     println("@test m.internalModel.linear_terms[lk$(i)][:y_type] == :($(m.internalModel.linear_terms[k][:y_type]))")
# end
#
# for i in 1:m.internalModel.num_constr_orig
#     if m.internalModel.constr_structure[i] == :affine
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:rhs] == $(m.internalModel.bounding_constr_mip[i][:rhs])")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:vars] == $(m.internalModel.bounding_constr_mip[i][:vars])")
#         println("@test isapprox(m.internalModel.bounding_constr_mip[$(i)][:coefs],$(m.internalModel.bounding_constr_mip[i][:coefs]);atol=1e-3)")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:sense] == :($(m.internalModel.bounding_constr_mip[i][:sense]))")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:cnt] == $(m.internalModel.bounding_constr_mip[i][:cnt])")
#     elseif m.internalModel.constr_structure[i] == :convex
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:rhs] == $(m.internalModel.bounding_constr_mip[i][:rhs])")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:vars] == $(m.internalModel.bounding_constr_mip[i][:vars])")
#         println("@test isapprox(m.internalModel.bounding_constr_mip[$(i)][:coefs],$(m.internalModel.bounding_constr_mip[i][:coefs]);atol=1e-3)")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:sense] == :($(m.internalModel.bounding_constr_mip[i][:sense]))")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:cnt] == $(m.internalModel.bounding_constr_mip[i][:cnt])")
#     else
#         error("Unknown constr_structure type $(m.constr_structure[i])")
#     end
# end
