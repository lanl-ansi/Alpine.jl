using Base.Test
using JuMP, Ipopt, MathProgBase
using POD

# Expression Testing Instances
include("../examples/bi1.jl")
include("../examples/exprs.jl")

# NLP Testing Instances
include("../examples/nlp1.jl")
include("../examples/nlp2.jl")
include("../examples/nlp3.jl")
include("../examples/castro02m2.jl")
include("../examples/castro03m2.jl")

# MINLP Testing Instances
include("../examples/util.jl")
include("../examples/meanvarx.jl")
include("../examples/blend029.jl")

include("solver.jl")
include("expression.jl")
include("algorithm.jl")
