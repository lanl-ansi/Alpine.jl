# Critical Updates for Alpine.jl
## v0.1.16
* Fixed CBC version to make it compatible 
* Removed support for "Mini" Formulation for convex hulls
* Updated tests and removed a few 

## v0.1.15
* Bug fixes for issues 134 and 135 (fractional exponents, constants within integer exponents)
* Added function 'expr_isolate_const' to nlexpr.jl
* Meaningful optimality gaps for zero upper bound values (issue 103)
* Typos in comments  
* Dropped support for colorful_alpine
* Ole Kroger's updates to docs/make.jl

## v0.1.10
* Multiple bug fixes in expression parsing 
* Added support for feasibility problems with constant objectives
* Typos in comments  

## v0.1.9
* Drop support for Julia < v1.0. 
* Move to the new Registrator system.
* The last supported version for Julia v0.6 is v0.1.8 
* The last supported version for Julia v0.7 is v0.1.7 (the bugs fixed in v0.1.8 is not ported to Julia v0.7)


## v0.1.8 
* Bug fixes in solver recognition 
* Updated Issue template and Readme.  
* Rounding bug fixed for binaries. 
* This version is only available for Julia v0.6

## v0.1.7
* Bug fixes for bounds on lifted variables representing nonlinear terms
* Inf bound warning issue fixed 

## v0.1.6
* Bug fix in operators.jl for accounting all extreme points in the bound evaluation 

## v0.1.5
* Unit exponent bug-fix
* Meaningful error messages for unsupported functions 
* Pinned tests to JuMP 0.18.5
* Cleaner & meaningful logs

## v0.1.3
* Finite large bounds on unbounded variables 

## Staged
* typo fix
* set Big-M values for unbounded variables

## v0.1.2
* Bug fix in bound tightening tolerance corner case

## v0.1.1
* Bug fix in Cplex solver recognition

## v0.1.0
* Alpine.jl first release

## Merged Features
* `disc_ratio` selection algorithm
* logarithm embedding formulation
* bounding model warm starting and no-good-cuts using solution pool
* `disc_vars` selection algorithm based on weighted minimum vertex cover
* more support for binary problem
* rounding heuristic and multi-start heuristic for feasible solutions
* recursive linear lifting
* some attempts to reason bounds for variables with infinite bounds

## New Operators
* `BINPROD`: product of binary vars
* `BINLIN`: product of a binary var and continuous var

## Data Structure
* Some refactoring to shorten the attribute names
* Eliminated unnecessary/obsolete user options
* Strengthen the checks to problems

## Code Structure
* `Alpine.jl`: "header" file
* `solver.jl`: data structure and problem loading
* `algorithm.jl`: high-level algorithms and major utilities (`presolve()`, `local_solve()`, `bounding_solve()`)
* `amp.jl`: algorithms for bounding formulation constructions
* `presolve.jl`: algorithms for presolves
* `heuristics.jl`: algorithms for heuristics
* `bounds.jl`: algorithms and utilities for handling bounds including bounds propagation and reasoning with lifting
* `mpb2moi.jl`: preparation for migration to `MOI` interface, a collection of all solver & modeling interface warpping functions
* `multi.jl`: high-level piecewise convex relaxation formulation construction and utilities for multilinear/monomial terms
* `tmc.jl`: old McCormick-based methods
* `const.jl`: (experimental) tries to categorize operators for cleaner and safer operations
* `embedding.jl`: logarithm embedding formulations and utilities
* `nlexpr.jl`: algorithms for expression handling
* `operators.jl`: operator detection algorithm and associated utilities for bounds, discretization variable selection, etc.
* `log.jl`: logging utilities
* `utility.jl`: general utility functions


## Debugs
* Fixed multiple bugs in recursive expression handling
* Fixed duplicates in constructing convex relaxations
* Fixed bound reasoning steps for complex terms
* Fixed many other things

## Tests & Coverage
* Most new features are covered (embedding, `disc_ratio` algorithm, etc.)
* Stronger expression parsing tests with new operators covered
* Coverage should be more than `~90%`
