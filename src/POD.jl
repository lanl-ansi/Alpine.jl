module POD

using JuMP
using MathProgBase
using Ipopt #Temporary Added

# package code goes here
include("expr.jl")
include("solver.jl")
include("algorithm.jl")

# Developing scripts
include("../dev/model.jl")

end # module
