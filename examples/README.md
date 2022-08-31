# Examples for Alpine
`run_examples.jl` file contains the main driver file which can be executed to compute global optimality for a non-convex MINLP of your choice. Multiple such JuMP models of MINLPs can be found inside the `examples/MINLPs` folder.

A few useful *hints* to explore within Alpine: 
* Try different integer values (>= 4) for `disc_ratio` to potentially observe 
    better Alpine run times, which can be instance specific.
* Choose `presolve_bt` to `false` if you prefer the optimization-based bound tightening (OBBT) presolve to be turned off. 
* If you prefer to use Alpine for only OBBT presolve, without any paritioning applied to the 
    nonlinear terms, include option `"apply_partitioning"` and set it to `false` while initializing the 
    Alpine solver. 