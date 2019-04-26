# Alpine, A global solver for nonconvex MINLPs <span style="color:black"></span>

Dev: [![Build Status](https://travis-ci.org/lanl-ansi/Alpine.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/Alpine.jl)
[![codecov](https://codecov.io/gh/lanl-ansi/Alpine.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/Alpine.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://lanl-ansi.github.io/Alpine.jl/latest/)

"ALPINE: glob(AL) o(P)timization for mixed-(I)nteger programs with (N)onlinear (E)quations", formerly **POD.jl**, is a novel global optimization solver that uses an adaptive, piecewise convexification scheme and constraints programming methods to solve Mixed-Integer Non-Linear Programs (non-convex MINLPs) efficiently. MINLPs are famously known as the "hard" programming problems that exist in many applications (see [MINLPLib.jl](https://github.com/lanl-ansi/MINLPLib.jl)). Alpine is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadradic Convex Programming (MIQCP), Non-Linear Programming (NLP), etc.

Unlike many other state-of-the-art MINLP solvers, Alpine is entirely built upon [JuMP](https://github.com/JuliaOpt/JuMP.jl) and [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) Interface in Julia, which provides incredible flexibility for usage and further development.

Alpine globally solves a given MINLP by:

* Analyzing the problem's expressions (objective & constraints) and applies approporite convex relaxations

* Performing novel adaptive partitioning methods to create piecewise relaxations, bound tightening and polyhedral outer-approximations to guarantee global convergence

**Allowable nonlinearities**: Alpine can currently handle MINLPs with polynomials in constraints and/or in the objective. 

**Alpine currently does not support problems with**: Exponential cones, Positive Semi-Definite (PSD) cones, Nonlinear functions including exp(x), ln(x), x/y, abs(x), x^y, a^x, trigonometric functions (sin(x), cos(x), etc), where 'a' is a constant, and 'x' and 'y' are variables. 

<!-- 
 **Illustration of Alpine's dynamic partitioning and outer-approximation on simple functions** ([Source](https://arxiv.org/abs/1707.02514))
 
<p align="center"> <img src="https://github.com/lanl-ansi/Alpine.jl/blob/master/Dynamic_partitions_github.png" width="580" class="centerImage"> </p>
-->
**Presentation on Alpine.jl at the [2nd Annual JuMP-dev Workshop](http://www.juliaopt.org/meetings/bordeaux2018/), held at the Institut de Mathématiques de Bordeaux, June 2018** 

[<img src="https://github.com/lanl-ansi/Alpine.jl/blob/master/alpine_slide.png" width="600" height="350">](https://www.youtube.com/watch?v=mwkhiEIS5JA)

## Installation

Alpine, with it's repository under the LANL-ANSI group, can be installed through the Julia package manager:

`julia> Pkg.add("Alpine")`

Developers: Any further development of Alpine can be conducted on a new branch or a forked repo.

## Underlying solvers

Though the algorithm implemented in Alpine is quite involved, most of the hard work and computational bottleneck would arise in the underlying solvers. Since every iteration of Alpine solves a subproblem to optimality, which is typically a convex MILP/MIQCQP and solves a nonconvex NLP/MINLP to local optimality, Alpine's run time heavily depends on the quality of these solvers. For best performance of Alpine, use commercial solvers such as CPLEX/Gurobi. However, due to the flexibility offered by JuMP/MathProgBase, the following solvers are supported in Alpine: 


| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)             |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Clp.jl)                 |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)           |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/JuliaOpt/AmplNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)           |

As the development of Alpine continues, support for solvers like [Mosek](http://www.mosek.com/), [Pajarito](https://github.com/JuliaOpt/Pajarito.jl), [GLPK](http://www.gnu.org/software/glpk/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), [Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite) is in the roadmap.

## Bug reports and support
Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

[issue tracker]: https://github.com/lanl-ansi/Alpine.jl/issues

## Challenging MINLPs
We are seeking out hard benchmark instances for MINLPs. Please get in touch either by opening an issue or [privately](https://harshangrjn.github.io/#contact) if you would like to share any hard instances. Challenging problems will help us determine how to improve Alpine.

## Citing Alpine

If you find Alpine useful in your work, we kindly request that you cite the following papers ([pdf](http://harshangrjn.github.io/pdf/JOGO_2018.pdf), [pdf](http://harshangrjn.github.io/pdf/CP_2016.pdf))
```
@article{NagarajanLuWangBentSundar2019,
  author = {Nagarajan, Harsha and Lu, Mowen and Wang, Site and Bent, Russell and Sundar, Kaarthik},
  title = {An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs},
  journal = {Journal of Global Optimization},
  year = {2019},
  issn = {1573-2916},
  doi = {10.1007/s10898-018-00734-1},
}

@inproceedings{NagarajanLuYamangilBent2016,
  title = {Tightening {McC}ormick relaxations for nonlinear programs via dynamic multivariate partitioning},
  author = {Nagarajan, Harsha and Lu, Mowen and Yamangil, Emre and Bent, Russell},
  booktitle = {International Conference on Principles and Practice of Constraint Programming},
  pages = {369--387},
  year = {2016},
  organization = {Springer},
  doi = {10.1007/978-3-319-44953-1_24},
}
```
