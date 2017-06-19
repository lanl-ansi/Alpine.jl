# POD, the solver

# Basic Usage with example of nlp1
`Pkg.clone("https://github.com/lanl-ansi/POD.git")`

`m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(),timeout=40, rel_gap=0.1, discretization_var_pick_algo=0)`

# NEED TO KNOW
1. The algorithm is currently under development.
2. Please limited usage of pranthese on constraints where non-linear terms exists.
3. Currently, the algorithm only handles simpler form of the bilinear formulation.

# Current Development
## Basic Solver Layout
Solver has a PODNonlinearModel structure that stores info utilized in the algorithm
## Basic POD Problem Type
TODO
## Simple Manipulation
TODO
## Local Solve
TODO

# Test Cases
## Expression Parsing
bi1, expr, nlp3,
## Solving Testing
nlp1, nlp3, util, castroEx02m2

See more in "examples" folder

# Unit Test
## Usage
After installation do `Pkg.test("POD")`
## Developed
Basic bilinear expression parsing/lifting/affine conversion test cases built
## TODO
1. Build some initialization testing for a clean POD model
2. Build more complex expressions for testing
3. Build basic solve-objective checking tests
4. Build formulation construction tests on tighten mccormicks
5. Build functional features tests
6. Build solver parameter tests

# More documentation to come
