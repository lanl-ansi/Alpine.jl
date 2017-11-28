# Usage

## Building a problem
For simplicity, we present a simple nonlinear program without any discrete variable for demonstration. The model is implemented in JuMP format. Two opon-source sub-solvers, [Ipopt](https://github.com/JuliaOpt/Ipopt.jl) and [Cbc](https://github.com/JuliaOpt/Cbc.jl), are used to handle local NLP (`nlp_local_solver`) and relaxed MIP (`mip_solver`). We control POD.jl to provide basic solving log through `log_level`. Before you solve, make sure that the two sub-solvers are properly installed and built.

```
using JuMP, POD
using Ipopt, Cbc

m = Model(solver=PODSolver(nlp_local_solver = IpoptSolver(print_level=0),
                           mip_solver = CbcSolver(logLevel=0),
                           rel_gap = 0.001,
                           log_level = 1))

@variable(m, x[1:8])

setlowerbound(x[1], 100)
setlowerbound(x[2], 1000)
setlowerbound(x[3], 1000)
setlowerbound(x[4], 10)
setlowerbound(x[5], 10)
setlowerbound(x[6], 10)
setlowerbound(x[7], 10)
setlowerbound(x[8], 10)

setupperbound(x[1], 10000)
setupperbound(x[2], 10000)
setupperbound(x[3], 10000)
setupperbound(x[4], 1000)
setupperbound(x[5], 1000)
setupperbound(x[6], 1000)
setupperbound(x[7], 1000)
setupperbound(x[8], 1000)

@constraint(m, 0.0025*(x[4]+x[6]) <= 1)
@constraint(m, 0.0025*(x[5]-x[4]+x[7]) <= 1)
@constraint(m, 0.01(x[8]-x[5]) <= 1)
@NLconstraint(m, 100*x[1] - x[1]*x[6] + 833.33252*x[4] <= 83333.333)
@NLconstraint(m, x[2]*x[4] - x[2]*x[7] - 1250*x[4] + 1250*x[5] <= 0)
@NLconstraint(m, x[3]*x[5] - x[3]*x[8] - 2500*x[5] + 1250000 <= 0)

@objective(m, Min, x[1]+x[2]+x[3])

solve(m)
```

## Interpreting POD log

POD.jl will output the log below. The output contains two main sections: 1) basic problem dimension and solver configuration summaries, and 2)presolve and main algorithms process. In the latter, The first two columns `NLP` and `MIP` are local solves and realxation solve objective value. Column `Objective` and `Bound` tracks the incumbent upper and lower bound of the problem. Without the loss of generality, upper bound is always defined through a feasible solution. Column `GAP%` is relative optimality gap evaluated through,

```math
    \textbf{Gap%} = \frac{|UB-LB|}{ϵ+|UB|} \times 100\%
```

Total wall time performance is tracked through column `CLOCK`. Last but not least, the iteration is tracked incrementally in column `Iter`.

```
full problem loaded into POD
problen sense Min
number of constraints = 6
number of non-linear constraints = 3
number of linear constraints = 3
number of variables = 8
MINLP local solver = POD
NLP local solver = Ipopt
MIP solver = Cbc
maximum solution time = Inf
maximum iterations =  99
relative optimality gap criteria = 0.00100 (0.1000 %)
detected nonlinear terms = 5
number of variables involved in nonlinear terms = 8
number of selected variables to discretize = 8
bilinear treatment = convex hull formulation
monomial treatment = convex hull formulation
using convex hull : sos2 formulation
using method 0 for picking discretization variable...
adaptively adding discretization ratio = 4

POD algorithm presolver started.
Performing local solve to obtain a feasible solution.
Presolve ended.
Presolve time = 0.02s
m.logs[:time_left] = Inf
 | NLP           | MIP           || Objective     | Bound         | GAP%          | CLOCK         | | Iter
 | 7049.2479     | 3004.247      || 7049.2479     | 3004.247      | 57.38202      | 0.07s         | | 1
 | -             | 4896.6075     || 7049.2479     | 4896.6075     | 30.53716      | 0.33s         | | 2
 | 7049.2479     | 5871.5307     || 7049.2479     | 5871.5307     | 16.70699      | 0.75s         | | 3
 | -             | 6717.2923     || 7049.2479     | 6717.2923     | 4.70909       | 2.8s          | | 4
 | 7049.267      | 6901.8131     || 7049.2479     | 6901.8131     | 2.0915        | 8.55s         | | 5
 | 7049.9758     | 7020.723      || 7049.2479     | 7020.723      | 0.40465       | 15.81s        | | 6
 | 7050.6454     | 7042.4094     || 7049.2479     | 7042.4094     | 0.09701       | 31.68s        | | 7
 | 7050.7924     | 7042.4094     || 7049.2479     | 7042.4094     | 0.09701       | 31.69s        | | finish

 POD ended with status Optimal
:Optimal
```
