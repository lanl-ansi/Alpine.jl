# Choosing Sub-Solvers

Solvers and solver-specific parameters are specified by POD, which are provided by particular solver packages. POD requires both a Non-linear Program Local Solver and a Mixed-integer Program Solver to work. For example, the `Clp` package exports a `ClpSolver` object, which can be passed to `PODSolver()` as follows::

```julia
    using JuMP
    using POD
    using Gurobi, Ipopt
    m = Model()
    setsolver(m, PODSolver(nlp_local_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(OutputFlag=0)))
```
