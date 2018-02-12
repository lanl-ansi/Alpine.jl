using Base.Test
using JuMP, MathProgBase
using Ipopt, Cbc, Pajarito
using GLPKMathProgInterface
using POD

poddir = Pkg.dir("POD")

# Expression Testing Instances
include("$(poddir)/test/examples/exprstest.jl")

# NLP Testing Instances
include("$(poddir)/test/examples/nlp.jl")
include("$(poddir)/test/examples/castro.jl")

# Multilinear Testing Instances
include("$(poddir)/test/examples/blend.jl")
include("$(poddir)/test/examples/multi.jl")

# Special operators
include("$(poddir)/test/examples/discretemulti.jl")
include("$(poddir)/test/examples/div.jl")
include("$(poddir)/test/examples/circle.jl")
include("$(poddir)/test/examples/convex.jl")
include("$(poddir)/test/examples/linearlift.jl")
include("$(poddir)/test/examples/specialopts.jl")
include("$(poddir)/test/examples/brainpc3.jl")

# Performe Tests
include("$(poddir)/test/solver.jl")
include("$(poddir)/test/expression.jl")
include("$(poddir)/test/algorithm.jl")
include("$(poddir)/test/utility.jl")

# Test Generation Code commented
# for k in keys(m.internalModel.nonlinear_terms)
#     println("@test m.internalModel.nonlinear_terms[$(k)][:y_idx] == $(m.internalModel.nonlinear_terms[k][:y_idx])")
#     println("@test m.internalModel.nonlinear_terms[$(k)][:id] == $(m.internalModel.nonlinear_terms[k][:id])")
#     println("@test m.internalModel.nonlinear_terms[$(k)][:lifted_constr_ref] == :($(m.internalModel.nonlinear_terms[k][:lifted_constr_ref]))")
#     println("@test m.internalModel.nonlinear_terms[$(k)][:nonlinear_type] == :$(m.internalModel.nonlinear_terms[k][:nonlinear_type])")
#     println("@test m.internalModel.nonlinear_terms[$(k)][:y_type] == :$(m.internalModel.nonlinear_terms[k][:y_type])")
#     println("@test m.internalModel.nonlinear_terms[$(k)][:constr_id] == $(m.internalModel.nonlinear_terms[k][:constr_id])")
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
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:coefs] == $(m.internalModel.bounding_constr_mip[i][:coefs])")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:sense] == :($(m.internalModel.bounding_constr_mip[i][:sense]))")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:cnt] == $(m.internalModel.bounding_constr_mip[i][:cnt])")
#     elseif m.internalModel.constr_structure[i] == :convex
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:rhs] == $(m.internalModel.bounding_constr_mip[i][:rhs])")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:vars] == $(m.internalModel.bounding_constr_mip[i][:vars])")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:coefs] == $(m.internalModel.bounding_constr_mip[i][:coefs])")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:sense] == :($(m.internalModel.bounding_constr_mip[i][:sense]))")
#         println("@test m.internalModel.bounding_constr_mip[$(i)][:cnt] == $(m.internalModel.bounding_constr_mip[i][:cnt])")
#     else
#         error("Unknown constr_structure type $(m.constr_structure[i])")
#     end
# end
