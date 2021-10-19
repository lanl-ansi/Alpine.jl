# Alpine's Documentation

```@meta
CurrentModule = Alpine
```
## Overview
[Alpine.jl](https://github.com/lanl-ansi/Alpine.jl) is a Julia/JuMP-based global optimization solver for noncovex programs. Alpine applies a two-stage approach to strengthen piecewise convex relaxations for mixed-integer nonlinear programs (MINLP) with multi-linear terms. In the first stage, Alpine exploits Constraint Programing techniques to contract the variable bounds. Further, it applies optimality-based bound tightening (OBBT) methods iteratively until a fixed point with respect to the bounds are achieved. In the second stage, Alpine partitions the variable domains using an adaptive multivariate partitioning scheme (a parametrized partition at the current solution of sequentially solved MILPs), leading to sparser but tighter relaxations. This approach decouples the number of partitions from the size of the variable domains, leading to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. Alpine also applies polyhedral cutting plane methods to handle convex relaxations of higher-order monomial terms.

**Allowable nonlinearities**: Alpine can currently handle MINLPs with polynomials in constraints and/or in the objective. Currently, there is no support for exponential cones and Positive Semi-Definite (PSD) cones in MINLPs. Alpine is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadratically Constrainted Quadradic Programs (MIQCQPs), Non-Linear Programs (NLPs), etc.

For more details, check out this [video](https://www.youtube.com/watch?v=mwkhiEIS5JA) on Alpine at the [2nd Annual JuMP-dev Workshop](http://www.juliaopt.org/meetings/bordeaux2018/), held at the Institut de Mathématiques de Bordeaux, June 2018.

## Installation
To use Alpine, first [download and install](https://julialang.org/downloads/) Julia. Note that the current version of Alpine is compatible with Julia 1.0 and later. 

The latest stable release of Alpine can be installed using the Julia package manager with

```julia
import Pkg
Pkg.add("Alpine")
```

At least one local nonlinear-programming (NLP) solver and a convex mixed-integer programming (MIP) solver are required for running Alpine. If your problem is an NLP, we recommend [Ipopt](https://github.com/jump-dev/Ipopt.jl), and else if your problem is a mixed-integer NLP (MINLP), we recommend [Juniper](https://github.com/lanl-ansi/Juniper.jl). For the underlying MIP solver, [Gurobi](https://github.com/jump-dev/Gurobi.jl) is highly recommended, as it is fast, scaleable and can be used to solve on fairly large-scale nonconvex programs. However, the open-source [CBC](https://github.com/jump-dev/Cbc.jl) or [GLPK](https://github.com/jump-dev/GLPK.jl) solver is also compatible with Alpine. For example, Ipopt and GLPK can be installed via the package manager with

```julia
import Pkg
Pkg.add("Ipopt")
Pkg.add("Gurobi")
```
For more details on choosing sub-solvers for Alpine, check out [here](https://lanl-ansi.github.io/Alpine.jl/latest/choosingsolver/). 

## Unit Tests
To run the tests in the package, run the following command after installing the Alpine package.

```julia
import Pkg
Pkg.test("Alpine")
```

## Citing Alpine
If you find Alpine useful in your work, we request you to cite the following paper [\[link\]](https://harshangrjn.github.io/pdf/JOGO_2018.pdf): 
```bibtex
@article{alpine_JOGO2019,
  author = {Nagarajan, Harsha and Lu, Mowen and Wang, Site and Bent, Russell and Sundar, Kaarthik},
  title = {An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs},
  journal = {Journal of Global Optimization},
  year = {2019},
  issn = {1573-2916},
  doi = {10.1007/s10898-018-00734-1},
}
```