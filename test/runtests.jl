using Logging
Logging.configure(level=ERROR)
using Base.Test
using JuMP, Ipopt, MathProgBase

include("expression.jl")
