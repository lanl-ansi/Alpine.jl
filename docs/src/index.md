# Alpine's Documentation

```@meta
CurrentModule = Alpine
```
## Overview
**[Alpine](https://github.com/lanl-ansi/Alpine.jl)** is a Julia/JuMP-based global optimization solver for non-convex programs. Alpine applies a two-stage approach to strengthen piecewise convex relaxations for mixed-integer nonlinear programs (MINLP) with multi-linear terms. In the first stage, Alpine exploits Constraint Programing techniques to contract the variable bounds. Further, it applies optimality-based bound tightening (OBBT) methods iteratively until a fixed point with respect to the bounds are achieved. In the second stage, Alpine partitions the variable domains using an adaptive multivariate partitioning scheme (a parametrized partition at the current solution of sequentially solved MILPs), leading to sparser but tighter relaxations. This approach decouples the number of partitions from the size of the variable domains, leading to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. Alpine also applies polyhedral cutting plane methods to handle convex relaxations of higher-order monomial terms.

**Allowable nonlinearities**: Alpine can currently handle MINLPs with polynomials in constraints and/or in the objective. Currently, there is no support for exponential cones and Positive Semi-Definite (PSD) cones in MINLPs. Alpine is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadratically Constrainted Quadradic Programs (MIQCQPs), Non-Linear Programs (NLPs), etc.

For more details, check out this [video](https://www.youtube.com/watch?v=mwkhiEIS5JA) on Alpine at the [2nd Annual JuMP-dev Workshop](http://www.juliaopt.org/meetings/bordeaux2018/), held at the Institut de Mathématiques de Bordeaux, June 2018.

## Installation
To use Alpine, first [download and install](https://julialang.org/downloads/) Julia. Note that the current version of Alpine is compatible with Julia 1.0 and later. 

The latest stable release of Alpine can be installed using the Julia package manager with

```julia
import Pkg
Pkg.add("Alpine")
```

At least one local nonlinear-programming (NLP) solver and a convex mixed-integer programming (MIP) solver are required for running Alpine. If your problem is an NLP, we recommend [Ipopt](https://github.com/jump-dev/Ipopt.jl), and else if your problem is a mixed-integer NLP (MINLP), we recommend [Juniper](https://github.com/lanl-ansi/Juniper.jl). For the underlying MIP solver, [Gurobi](https://github.com/jump-dev/Gurobi.jl) is highly recommended, as it is fast, scaleable and can be used to solve on fairly large-scale non-convex programs. However, the open-source [HiGHS](https://github.com/jump-dev/HiGHS.jl) or [GLPK](https://github.com/jump-dev/GLPK.jl) solver is also compatible with Alpine. For example, Ipopt and Gurobi can be installed via the package manager with

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
If you find Alpine useful in your work, we request you to cite the following paper [\[pdf\]](https://harshangrjn.github.io/pdf/JOGO_2018.pdf) [\[pdf\]](http://harshangrjn.github.io/pdf/CP_2016.pdf): 
```bibtex
@article{alpine_JOGO2019,
  title = {An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs},
  author = {Nagarajan, Harsha and Lu, Mowen and Wang, Site and Bent, Russell and Sundar, Kaarthik},
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
If you find the underlying piecewise polyhedral formulations implemented in Alpine useful in your work, we kindly request that you cite the following papers ([link-1](https://doi.org/10.1016/j.orl.2020.12.002), [link-2](http://www.optimization-online.org/DB_HTML/2022/07/8974.html)): 
```bibtex
@article{alpine_ORL2021,
  title = {Piecewise polyhedral formulations for a multilinear term},
  author = {Sundar, Kaarthik and Nagarajan, Harsha and Linderoth, Jeff and Wang, Site and Bent, Russell},
  journal = {Operations Research Letters},
  volume = {49},
  number = {1},
  pages = {144--149},
  year = {2021},
  publisher = {Elsevier}
}

@article{alpine_OptOnline2022,
  title={Piecewise Polyhedral Relaxations of Multilinear Optimization},
  author={Kim, Jongeun and Richard, Jean-Philippe P. and Tawarmalani, Mohit},
  eprinttype={Optimization Online},
  date={2022}
}
```

If you find the random families of non-convex QCQP instances (using generator files within this [folder](https://github.com/lanl-ansi/Alpine.jl/tree/master/examples/random_QCQPs)) useful for your research, we kindly request that you cite the following paper [arXiv link](https://arxiv.org/abs/2301.00306):
```bibtex
@article{alpine_learning_2022,
  title={Learning to Accelerate the Global Optimization of Quadratically-Constrained Quadratic Programs},
  author={Kannan, Rohit and Nagarajan, Harsha and Deka, Deepjyoti},
  journal={arXiv preprint:2301.00306},
  url={https://arxiv.org/abs/2301.00306},
  year={2022}
}
```