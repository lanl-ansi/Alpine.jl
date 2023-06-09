# Examples for Alpine
`run_examples.jl` file contains the main driver file which can be executed to compute global optimality for a non-convex MINLP of your choice. Multiple such JuMP models of MINLPs can be found inside the `examples/MINLPs` folder.

A few useful **hints** to explore within Alpine: 
* Try different integer values (>= 4) for `partition_scaling_factor` to potentially observe 
    better Alpine run times, which can be instance specific.
* Choose `presolve_bt` to `false` if you prefer the optimization-based bound tightening (OBBT) presolve to be turned off. 
* If you prefer to use Alpine for only OBBT presolve, without any paritioning applied to the 
    nonlinear terms, include option `"apply_partitioning"` and set it to `false` while initializing the 
    Alpine solver. 
* If the nonconvex instance you want to solve has binary variables (i.e., a MINLP like in [here](https://github.com/lanl-ansi/Alpine.jl/blob/32757b2568acfe71257e41ed4122a120717a6460/examples/MINLPs/blend.jl#LL1C1-L1C1), then include a local MINLP solver as an additional attribute while initializing `Alpine.Optimizer` in `run_examples.jl`. 

# Random families of non-convex QCQPs
Scripts for generating homogeneous families of nonconvex QCQP instances, useful for benchmarking ML-based models for tuning certain algorithmic parameters of global optimization, can be found [here](https://github.com/lanl-ansi/Alpine.jl/tree/master/examples/random_QCQPs). 