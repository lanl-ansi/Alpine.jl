using JuMP, MathOptInterface
using Ipopt, Cbc
using GLPK
using Alpine

using Test

using Random
using Pkg
alpine_dir = joinpath(dirname(pathof(Alpine)), "..")

@test true == true