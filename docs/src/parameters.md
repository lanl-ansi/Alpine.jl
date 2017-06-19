# Parameters

```@meta
CurrentModule = POD
```

## General

There are some general paremters that controls the behavior of the AMP:

* `log_level`: verbosity of the algorihtm (`default=0`), set to `1` for turning on logging, `100` for detailed debugging mode
* `timeout`: total time regulated for the solver (`default=Inf`), in seconds
* `maxiter`: total iteration allowed in the [`global_solve`](@ref) (`default=999`)
* `rel_gap`: relative gap considered for optimality during [`global_solve`](@ref) (`default=1e-4`), evaluated by ``\frac{UB-LB}{UB} 100 \times \%``
* `tolerance`: numerical tolerance used furing the process of [`global_solve`](@ref) (`default=1e-6`)

## Discretization-based parameters

More parameter descriptions to come...
