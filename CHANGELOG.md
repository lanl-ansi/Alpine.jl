# Critical Updates for Alpine.jl

I've forked the master branch hosted by lanl-ansi to my own Github while keeping all the most recent features merged on my master. Current, I am testing MPC instances with my master branch. Here is a summary of what are changed:

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
