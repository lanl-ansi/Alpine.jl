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
minmax_bound_tightening
create_bound_tightening_model
solve_bound_tightening_model
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
