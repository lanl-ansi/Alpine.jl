# Parameters

## General Solver Control
| Keyword    | Type | Default | Description |
|------------|------|---------|-------------|
| `log_level`   | `Int`     | `1`       | Option to control log output. `0`: not log, `1`: standard log, `100`: debug log. |
| `timeout`     | `Float`   | `Inf`     | Maximum wall time allowed in seconds before termination. |
| `max_iter`    | `Int`     | `99`      | Maximum iteration of the main algorithm allowed before termination. |
| `rel_gap`     | `Float`   | `0.0001`  | Optimal relative gap tolerance before termination. |
| `tol`         | `Float`   | `1e-6`    | General numerical tolerance used to resolve numerical issues. |
| `fea_tol`     | `Float`   | `1e-8`    | Feasibility tolearnace used to resolve numerical issues.

## Sub-solvers Control
| Keyword               | Type    | Default   | Description |
|-----------------------|---------|-----------|-------------|
| `minlp_local_solver`  | `AbstractMathProgSolver` | `UnsetSolver` | MINL sub-solver to solve any MINLP locally. This solver is not required when solving NLP. |
| `nlp_local_solver`    | `AbstractMathProgSolver` | `UnsetSolver` | NLP sub-solver to solve any NLP locally. This solver is required. It missing, it will use `minlp_local_solver`'s options. |
| `mip_solver`          | `AbstractMathProgSolver` | `UnsetSolver` | MILP sub-solver to solve any MILP/IP/LP/MIQCP globally. This solver is required. |

## Main Algorithm Control
| Keyword               | Type      | Default | Description |
|-----------------------|-----------|---------|-------------|
| `presolve`            | `Bool`    | `true`  | Option to turn on and off the presolve algorithm for bound tightening |
| `presolve_algo`       | `Int`     | `1`     | Algorithm used in presolve. `0`: same as turn presolve off, `1`: simple optimality-based bound tightening algorithm using root relaxation, `2`: optimality-based bound tightening algorithm using tighter relaxation. |
| `presolve_max_iter`   | `Int`     | `99`    | Maximum ieration allowed for any sequential bound tightening algorithm.
| `presolve_timeout`    | `Int`     | `Inf`   | Maximum time allowed for each relaxation optimization solved during bound tightening algorithm. |
| `recognize_convex`    | `Bool`    | `false` | (Beta) Option to ask the solver to recognize built-in patterns of convex constraints for stronger relaxation. |
| `disc_ratio`          | `Int`     | `8`     | Discretization ratios used when dynamically building variable partitions. Value must be `>2`. |
| `disc_var_pick_algo`  | `Int`     | `2`     | Option to choose which variables' domain should be partitioned in the solver. `0`: choose all variables in nonlinear terms, `1`: use a minimum vertex cover to choose the smallest variable set for algorithm convergence, `2`: balace between `0` and `1` given the problem size, `3`: heuristic algorithm to dynammically update the choice during algorithmic process. |
| `disc_abs_width_tol`  | `Float`   | `1e-4` | Option to control the minimum absolute length of a partition in the variable domain. No further partition will be added to very small partitioned domains for global convergence. |
| `disc_rel_width_tol`  | `Float`   | `1e-6` | Option to control the minimum relative length of a partition to the original variable domain. No further partition will be added to very small partitioned domains for global convergence.  |
| `disc_consec_forbid`  | `Int`     | `3`    | Option to control maximum encoutner of a same relaxation solution consecutively to avoid local trap. New partitions will be added elsewhere if this is triggered. Use `0` to turn this feature off. |
| `bounds_propagation`  | `Bool`    | `true` | Options to control the basic bound propagation to reason each variable's bounds. POD is currently with limited support for unbounded problem. |

## User API Control (Advanced Control)
| Keyword                   | Type      | Default   | Description |
|---------------------------|-----------|-----------|-------------|
| `nlterm_patterns`         | `Any`     | `nothing` | Methods to recognize user-defined patterns for nonlinear term during expression parsing. User-inputs override built-in methods. |
| `convex_patterns`         | `Any`     | `nothing` | Methods to recognize user-defined convex constraints patterns when `recognize_convex` is true. User-inputs override built-in methods. |
| `convexify_methods`       | `Any`     | `nothing` | Methods to convexify nonlinear terms during construction of the convex relxation. User-inputs override built-in methods. |
| `disc_partition_method`   | `Any`     | `nothing` | Methods to add partitions during the main algorithm process. User-inputs override built-in methods. Change of this default option may affect the global convergence. |
