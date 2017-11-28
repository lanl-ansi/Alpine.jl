# Parameters

## General Solver Control
| Keyword    | Type | Description |
|---------------|-----------|-----|
| `log_level`   | `Int`     | ... |
| `timeout`     | `Float`   | ... |
| `max_iter`    | `Int`     | ... |
| `rel_gap`     | `Float`   | ... |
| `tol`         | `Float`   | ... |

## Sub-solvers Control
| Keyword               | Type             | Description |
|-----------------------|--------------------------|-----|
| `minlp_local_solver`  | `AbstractMathProgSolver` | ... |
| `nlp_local_solver`    | `AbstractMathProgSolver` | ... |
| `mip_solver`          | `AbstractMathProgSolver` | ... |

## Main Algorithm Control
| Keyword            | Type | Description |
|-----------------------|-----------|-----|
| `presolve`            | `Int`     | ... |
| `presolve_algo`       | `Int`     | ... |
| `presolve_max_iter`   | `Int`     | ... |
| `presolve_mip_relax`  | `Int`     | ... |
| `presolve_timeout`    | `Int`     | ... |
| `recognize_convex`    | `Bool`    | ... |
| `disc_ratio`          | `Int`     | ... |
| `disc_var_pick_algo`  | `Int`     | ... |
| `disc_abs_width_tol`  | `Int`     | ... |
| `disc_consec_forbid`  | `Int`     | ... |
| `bounds_propagation`  | `Bool`    | ... |

## User API Control
| Keyword                | Type | Description |
|---------------------------|-----------|-----|
| `nlterm_patterns`         | `Any`     | ... |
| `convex_patterns`         | `Any`     | ... |
| `convexify_methods`       | `Any`     | ... |
| `disc_partition_method`   | `Any`     | ... |
