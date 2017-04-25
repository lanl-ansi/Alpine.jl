module POD

using JuMP
using MathProgBase

# package code goes here
include("types.jl")
include("operators.jl")
include("expr.jl")
include("solver.jl")
include("algorithm.jl")

end # module
