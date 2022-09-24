# Functions

```@meta
CurrentModule = Alpine
```

## High-level Algorithmic Functions
These are the high-level algorithmic functions:
```@docs
presolve
global_solve
local_solve
bounding_solve
```

## Adapative Partitioning Methods
```@docs
create_bounding_mip
pick_disc_vars
fix_domains
min_vertex_cover
```

## Presolve Methods
```@docs
bound_tightening
optimization_based_bound_tightening
create_obbt_model
solve_obbt_model
resolve_var_bounds
```

## Utility Methods
```@docs
update_var_bounds
discretization_to_bounds
init_disc
to_discretization
flatten_discretization
add_adaptive_partition
```
