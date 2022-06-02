# Alpine, A global solver for non-convex MINLPs <span style="color:black"></span>

STATUS: [![CI](https://github.com/lanl-ansi/Alpine.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/lanl-ansi/Alpine.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/lanl-ansi/Alpine.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/Alpine.jl)
[![Documentation](https://github.com/lanl-ansi/Alpine.jl/actions/workflows/documentation.yml/badge.svg)](https://lanl-ansi.github.io/Alpine.jl/latest/)

"ALPINE: glob(AL) o(P)timization for mixed-(I)nteger programs with (N)onlinear (E)quations", is a novel global optimization solver that uses an adaptive, piecewise convexification scheme and constraint programming methods to solve non-convex Mixed-Integer Non-Linear Programs (MINLPs) efficiently. MINLPs are typically "hard" optimization problems which appear in numerous applications (see [MINLPLib.jl](https://github.com/lanl-ansi/MINLPLib.jl)). 

Alpine is entirely built upon [JuMP](https://github.com/jump-dev/JuMP.jl) and [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) in Julia, which provides incredible flexibility for usage and further development.

Alpine globally solves a given MINLP by:

* Analyzing the problem's expressions (objective & constraints) and applies appropriate convex relaxations and polyhedral outer-approximations

* Performing sequential optimization-based bound tightening (OBBT) and an iterative MIP-based adaptive partitioning scheme via piecewise polyhedral relaxations with a guarantee of global convergence

**Allowable nonlinearities**: Alpine can currently handle MINLPs with polynomials in constraints and/or in the objective. Currently, there is no support for exponential cones and Positive Semi-Definite (PSD) cones in MINLPs. Alpine is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadratically Constrainted Quadradic Programs (MIQCQPs), Non-Linear Programs (NLPs), etc.

<!-- 
 **Illustration of Alpine's dynamic partitioning and outer-approximation on simple functions** ([Source](https://arxiv.org/abs/1707.02514))
 
<p align="center"> <img src="https://github.com/lanl-ansi/Alpine.jl/blob/master/Dynamic_partitions_github.png" width="580" class="centerImage"> </p>
-->
For more details, check out this [video](https://www.youtube.com/watch?v=mwkhiEIS5JA) on Alpine.jl at the [2nd Annual JuMP-dev Workshop](http://www.juliaopt.org/meetings/bordeaux2018/), held at the Institut de Mathématiques de Bordeaux, June 2018.

<!-- [<img src="https://github.com/lanl-ansi/Alpine.jl/blob/master/alpine_slide.png" width="600" height="350">](https://www.youtube.com/watch?v=mwkhiEIS5JA) -->

## Installation and Usage

Alpine can be installed through the Julia package manager:

`julia> Pkg.add("Alpine")`

Developers: Any further development of Alpine can be conducted on a new branch or a forked repo.

Check the "test/examples" folder on how to use this package. 

## Underlying solvers

Though the MIP-based bounding algorithm implemented in Alpine is quite involved, most of the computational bottleneck arises in the underlying MIP solvers. Since every iteration of Alpine solves an MIP sub-problem, which is typically a convex MILP/MIQCQP, Alpine's run time heavily depends on the run-time of these solvers. For the best performance of Alpine, we recommend using the commercial solver [Gurobi](https://www.gurobi.com), which is avaible [free](https://www.gurobi.com/academia/academic-program-and-licenses/) for academic purposes. However, due to the flexibility offered by [JuMP](https://github.com/jump-dev/JuMP.jl), the following MIP and NLP solvers are supported in Alpine: 


| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)             |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/jump-dev/Cbc.jl)                 |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)           |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/jump-dev/AmplNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/jump-dev/KNITRO.jl)           |
| [Xpress](https://www.fico.com/en/products/fico-xpress-optimization)            | [Xpress.jl](https://github.com/jump-dev/Xpress.jl)
| [HiGHS](https://highs.dev/)            | [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl)

## Bug reports and support
Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

[issue tracker]: https://github.com/lanl-ansi/Alpine.jl/issues

## Challenging Problems
We are seeking out hard benchmark instances for MINLPs. Please get in touch either by opening an issue or [privately](https://harshangrjn.github.io/#contact) if you would like to share any hard instances.

## Citing Alpine

If you find Alpine useful in your work, we kindly request that you cite the following papers ([pdf](http://harshangrjn.github.io/pdf/JOGO_2018.pdf), [pdf](http://harshangrjn.github.io/pdf/CP_2016.pdf))
```bibtex
@article{alpine_JOGO2019,
  author = {Nagarajan, Harsha and Lu, Mowen and Wang, Site and Bent, Russell and Sundar, Kaarthik},
  title = {An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs},
  journal = {Journal of Global Optimization},
  year = {2019},
  issn = {1573-2916},
  doi = {10.1007/s10898-018-00734-1},
}

@inproceedings{alpine_CP2016,
  title = {Tightening {McCormick} relaxations for nonlinear programs via dynamic multivariate partitioning},
  author = {Nagarajan, Harsha and Lu, Mowen and Yamangil, Emre and Bent, Russell},
  booktitle = {International Conference on Principles and Practice of Constraint Programming},
  pages = {369--387},
  year = {2016},
  organization = {Springer},
  doi = {10.1007/978-3-319-44953-1_24},
}
```
