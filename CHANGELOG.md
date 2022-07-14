# Alpine.jl Change Log

## v0.4.1
- Added MOI attribute for tightened bounds post-OBBT (@blegat)

## v0.4.0
- Fixed OBBT convergence issue due to average bound reduction
- Added new termination criteria to OBBT: max_reduction, max_width, rel_gap at every iteration of OBBT
- Fixed the issue to tighten only original nonlinear variables (aux vars were counted before)
- `presolve_bt_bound_tol` changed to `1E-4` (numerically stable)
- `presolve_bt_obj_bound_tol` added to avoid infeasibility in OBBT iterations when bounds are tights
- Parameter tuning (tolerances) for stable OBBT convergence
- Fixed bug in `best_rel_gap` evaluation
- Bug fixed to check any user-defined operators, and disable hessians based on availabe features. Unit test added for this 
- Now, Alpine terminates (w/o invoking `global_solve`) when presolve finds the global optimum due to `rel_gap` eval at every OBBT iteration. Unit test added for this.
- Post-OBBT optimality gap (`presolve_best_rel_gap`) user-option added to be able to extract best opt gap after OBBT when Alpine terminates
- `apply_partitioning` feature added to be able to terminate Alpine without applying the MIP-based paritioning algorithm (mainly to run only OBBT)
- Revised unit test (in `test_algorithm`) for reduced run times

## v0.3.0
- JuMP v1.0+ in examples
- Pavito v0.3.5 support
- Dropped support for redundant `run_bounding_iteration` function 
- Clean-up and re-org of optimizer structs and `MOI_wrapper` functions
- Bug fix in `relax_integrality` under presolve
- Added Rosenbrock function
- Added user-definable `presolve_bt_improv_tol` for variable range reductions in OBBT
- Minor update on finite values for obj-bound in `create_bound_tightening_model`
- NLP solution rounding issue resolved (now compatible with default tol, 1e-6), to avoid the variable solution outside discretization (#190)     
- `circleN` instance updated to `circle_MINLPLib` to eliminate Pavito's mixed-integer cycling in unit tests (#190)
- Clean-up in printing variable solutions (`variable_values`)

## v0.2.10
- Add support for the HiGHS solver

## v0.2.9
- Update to JuMP v1.0

## v0.2.8
- Update to JuMP v0.22, v0.23
- Update to MOI v1.0
- Simplify test dependencies

## v0.2.7
- Support for module shortcut `Alp` in function calls
- Solver options clean up in `examples` folder
- Dropped support for unsued functions:
     `_bound_sense`
     `update_boundstop_options`
     `amp_post_inequalities_int`
     `amp_pick_ratevec`
     `amp_collect_tight_regions`
- Improved code coverage
- Docs cleanup

## v0.2.6
- Bug fix in `algorithm.jl` to handle Gurobi as the MIP solver with updated MOI 
- Deactivated redundant functions `collect_lb_pool`, `merge_solution_pool` and `track_new_partition_idx`
- Moved `examples` out of `test` folder
- Updated `run_examples` to handle Gurobi as the MIP solver
- Removed dependency on the `Distributed` package
- Improved code coverage
- Typos fixed in docs

## v0.2.5
- Cleaning up of solver logs
- Time correction for local solve 
- Removed unused `Alpine.Optimizer` options in solver.jl

## v0.2.4
- Added unit tests for utility functions
- Removed un-used functions from utility functions which aren't part of Alpine's algorithm     
- Dropped support for heuristic rounding algorithm (`heu_basic_rounding`) - `minlp_solver` input is a must
- Dropped support for adjusting branching priority in MIP solvers (`adjust_branch_priority`)


## v0.2.3
- Migrating from Travis to Github-actions
- Major documentation clean-up
- Dependency packages clean-up 
- runtests.jl clean-up
- Solver options update: `presolve_max_iter` => `presolve_bt_max_iter` and others
- Added `test/examples/run_example.jl`

## v0.2.2
- README clean up 

Closed issues: 
- Constraints are not supported and cannot be bridged into supported constrained variables and constraints (#164)

## v0.2.1
- Logging clean up 
- Removal of redundant tests

Closed issues: 
- Juniper's ALMOST_LOCALLY SOLVED error (#157)
- unncessary eval usages (#153)
- Logging issue with constants in objective (#159)

## v0.2.0
- MOI wrapper and JuMP v0.21 support 
- expression.jl clean up
- eval() removal where unnecessary (#153)
- Drop support for all trigonometric related functions (both in src and testing)
- Drop support of expression parsing for generic integer variables in MINLPs (both in src and testing)

Closed issues:
- arc consistency partitioning pruning scheme as an example extension for POD (#37)
- Generate convex relaxations for nonconvex MINLPs (#86)
- Error in convexification of monomials (#104)
- On variables chosen for partitioning (#121)
- Travis and Documeter Maintenance (#144)
- Compilation error (#146)
- bound propagation buggy (#147)
- Reenable test: Validation Test || PBT-AMP-TMC || basic solve || examples/nlp3.jl (#150)
- Reenable test: Validation Test || AMP-CONV-FACET || basic solve || examples/nlp1.jl (#151)

## v0.1.16
- Fixed CBC version to make it compatible 
- Removed support for "Mini" Formulation for convex hulls
- Updated tests and removed a few 

## v0.1.15
- Bug fixes for issues 134 and 135 (fractional exponents, constants within integer exponents)
- Added function 'expr_isolate_const' to nlexpr.jl
- Meaningful optimality gaps for zero upper bound values (issue 103)
- Typos in comments  
- Dropped support for colorful_alpine
- Ole Kroger's updates to docs/make.jl

## v0.1.10
- Multiple bug fixes in expression parsing 
- Added support for feasibility problems with constant objectives
- Typos in comments  

## v0.1.9
- Drop support for Julia < v1.0. 
- Move to the new Registrator system.
- The last supported version for Julia v0.6 is v0.1.8 
- The last supported version for Julia v0.7 is v0.1.7 (the bugs fixed in v0.1.8 is not ported to Julia v0.7)


## v0.1.8 
- Bug fixes in solver recognition 
- Updated Issue template and Readme.  
- Rounding bug fixed for binaries. 
- This version is only available for Julia v0.6

## v0.1.7
- Bug fixes for bounds on lifted variables representing nonlinear terms
- Inf bound warning issue fixed 

## v0.1.6
- Bug fix in operators.jl for accounting all extreme points in the bound evaluation 

## v0.1.5
- Unit exponent bug-fix
- Meaningful error messages for unsupported functions 
- Pinned tests to JuMP 0.18.5
- Cleaner & meaningful logs

## v0.1.3
- Finite large bounds on unbounded variables 

## Staged
- typo fix
- set Big-M values for unbounded variables

## v0.1.2
- Bug fix in bound tightening tolerance corner case

## v0.1.1
- Bug fix in Cplex solver recognition

## v0.1.0
- Alpine.jl first release

## Merged Features
- `disc_ratio` selection algorithm
- logarithm embedding formulation
- bounding model warm starting and no-good-cuts using solution pool
- `disc_vars` selection algorithm based on weighted minimum vertex cover
- more support for binary problem
- rounding heuristic and multi-start heuristic for feasible solutions
- recursive linear lifting
- some attempts to reason bounds for variables with infinite bounds

## New Operators
- `BINPROD`: product of binary vars
- `BINLIN`: product of a binary var and continuous var

## Data Structure
- Some refactoring to shorten the attribute names
- Eliminated unnecessary/obsolete user options
- Strengthen the checks to problems

## Code Structure
- `Alpine.jl`: "header" file
- `solver.jl`: data structure and problem loading
- `algorithm.jl`: high-level algorithms and major utilities (`presolve()`, `local_solve()`, `bounding_solve()`)
- `amp.jl`: algorithms for bounding formulation constructions
- `presolve.jl`: algorithms for presolves
- `heuristics.jl`: algorithms for heuristics
- `bounds.jl`: algorithms and utilities for handling bounds including bounds propagation and reasoning with lifting
- `mpb2moi.jl`: preparation for migration to `MOI` interface, a collection of all solver & modeling interface warpping functions
- `multi.jl`: high-level piecewise convex relaxation formulation construction and utilities for multilinear/monomial terms
- `tmc.jl`: old McCormick-based methods
- `const.jl`: (experimental) tries to categorize operators for cleaner and safer operations
- `embedding.jl`: logarithm embedding formulations and utilities
- `nlexpr.jl`: algorithms for expression handling
- `operators.jl`: operator detection algorithm and associated utilities for bounds, discretization variable selection, etc.
- `log.jl`: logging utilities
- `utility.jl`: general utility functions


## Debugs
- Fixed multiple bugs in recursive expression handling
- Fixed duplicates in constructing convex relaxations
- Fixed bound reasoning steps for complex terms
- Fixed many other things

## Tests & Coverage
- Most new features are covered (embedding, `disc_ratio` algorithm, etc.)
- Stronger expression parsing tests with new operators covered
- Coverage should be more than `~90%`
