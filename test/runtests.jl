using Logging
Logging.configure(level=ERROR)
using Base.Test
using JuMP, Ipopt, MathProgBase
using POD

include("solver.jl")
include("expression.jl")
