# POD, the solver

# Basic Usage
`Pkg.clone("https://github.com/lanl-ansi/POD.git")`

`using POD`

`m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0))`

# NEED TO KNOW
1. To be able to do `Pkg.test("POD")`, disable line 242-244 in JuMP->parsenlp.jl.

2. If not intend to test expression manipulation, it is not necessary.

# Current Development
## Basic Solver Layout
TODO
## Basic POD Problem Type
TODO
## Simple Manipulation
TODO
## Local Solve
TODO

# Test Cases
bi1, nlp1, nlp2, nlp3

See examples folder

# Test Usage
After installation do `Pkg.test("POD")`

# Limitations
1. Currently only handles bilinear constraints with limited usage of ()
2. Building model using @addNLconstraint() only

# More documentation to come
