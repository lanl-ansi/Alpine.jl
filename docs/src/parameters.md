# Parameters

```@meta
CurrentModule = Alpine
```

## General

There are some general parameters that controls the behavior of the AMP:

* `loglevel(default=0)`: verbosity of the algorithm, set to `1` for turning on logging, `100` for detailed debugging mode

* `timeout(default=Inf)`: total time regulated for the solver in seconds

* `maxiter(default=999)`: total iteration allowed in the [`global_solve`](@ref)

* `relgap(default=1e-4)`: relative gap considered for optimality during [`global_solve`](@ref). Bounds are evaluated using ``\frac{UB-LB}{UB} 100 \times \%``

* `tol(default=1e-6)`: numerical tol used during the process of [`global_solve`](@ref)

## Discretization-based parameters

* `disc_ratio(default=4)`: used during [`add_adpative_partition`](@ref) for measuring the radius of new partitioning relative to the active domain

* `disc_var_pick(default=0)`: controls algorithm/methods used for selecting variables for discretization, `0` is for max-cover, `1` is for minimum-vertex-cover. This parameter allows functional inputs.

* `disc_add_partition_method`: allows functional input on how new partitions should be constructed. This parameter is not fully stable.

## Presolve Parameters

* `presolve_track_time(default=false)`: consider presolve time as the total time or not

* `presolve_bt(default=false)`: perform built-in bound tightening presolve procedure

* `presolve_maxiter(default=9999)`: maximum iteration allowed using presolve process

* `presolve_bt_width_tol(default=1e-3)`: independent numerical tol used in presolve for bound tightening procedure. Note that this procedure is more sensitive to the tol in here. Small tol is more likely to results in strange presolve behavior.

* `presolve_bt_algo(default=1)`: method used to do built-in bound tightening, choose `1` for regular bounding tightening,  `2` for Tighten McCormick bound tightening.

* `presolve_bt_relax(default=false)`: whether to relax the bounding tightening MILP solved or not

* `presolve_bt_mip_timeout(default=Inf)`: time limit used for individual MILP solved during bound tightening presolve

More parameter descriptions to come...
