# Solver options for Alpine

```@meta
CurrentModule = Alpine
```

## General Options

Here are a few general solver options which control the performance of Alpine:

* `log_level(default=0)`: verbosity level of Alpine; choose `1` for turning on logging, else `100` for detailed debugging mode.

* `time_limit(default=Inf)`: total limit on run time for Alpine in seconds.

* `max_iter(default=999)`: total number of iterations allowed in the [`global_solve`](@ref). 

* `rel_gap(default=1e-4)`: relative gap considered for global convergence during [`global_solve`](@ref). Bounds are evaluated using ``\frac{UB-LB}{UB} \cdot 100 \%``.

* `tol(default=1e-6)`: numerical tolerance used during the process of [`global_solve`](@ref). 

## Adaptive Partitioning Options

* `disc_ratio(default=4)`: used during [`add_adaptive_partition`](@ref) for measuring the width of new partitions relative to the active partition chosen in the sequentially solved lower-bounding MIPs. This value can substantially affect the run time for global convergence; this value can be set to different integer values (>= 4) for various classes of problems. 

* `disc_var_pick(default=0)`: controls Alpine's algorithm used for selecting variables for partitioning; `0` is for max-cover, `1` is for minimum-vertex-cover. This parameter allows functional inputs.

* `disc_add_partition_method`: allows functional input on how new partitions could be constructed.

## Presolve Options

* `presolve_track_time(default=true)`: includes/excludes presolve run time in the total Alpine's run time. 

* `presolve_bt(default=false)`: performs sequential, optimization-based bound tightening (OBBT) at the presolve step. 

* `presolve_bt_max_iter(default=9999)`: maximum number of iterations allowed using the sequential OBBT step. 

* `presolve_bt_width_tol(default=1e-3)`: numerical tolerance value used in the OBBT step. Note that smaller values of this tolerance can lead to inaccurate solutions.  

* `presolve_bt_algo(default=1)`: method chosen to perform the presolve step; choose `1` for the built-in OBBT, else `2` for user-input functions for presolve.

* `presolve_bt_relax_integrality(default=false)`: relaxes the integrality of the existing integer variables in the OBBT step, if the input problem is an MINLP. 

* `presolve_bt_mip_time_limit(default=Inf)`: time limit for individual MILPs solved during the sequential OBBT procedure. 

Note that the above-mentioned list of solver options is not comprehensive, but can be found in [solver.jl](https://github.com/lanl-ansi/Alpine.jl/blob/master/src/solver.jl). 
