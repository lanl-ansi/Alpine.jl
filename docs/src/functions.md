# Functions

```@meta
CurrentModule = POD
```

## High-level Algorithmic Operations
These are the high-level algorithmic methods:
```@docs
presolve
global_solve
local_solve
bounding_solve
```

## Adapative Partitioning Methods
```@docs
create_bounding_mip
pick_vars_discretization
fix_domains
min_vertex_cover
max_cover
```

## Presolve Methods
```@docs
presolve_bounds_tightening
minmax_bounds_tightening
create_bound_tightening_model
solve_bound_tightening_model
resolve_lifted_var_bounds
```

## Utility Methods
```@docs
update_var_bounds
discretization_to_bounds
initialize_discretization
to_discretization
flatten_discretization
add_discretization
update_mip_time_limit
fetch_timeleft_symbol
```
