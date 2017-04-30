using Base.Test
using JuMP, Ipopt, MathProgBase
using POD

include("../examples/bi1.jl")
include("../examples/nlp1.jl")
include("../examples/nlp2.jl")
include("../examples/nlp3.jl")

include("solver.jl")
include("expression.jl")
