# Alpine, a global solver for non-convex MINLPs

[![CI](https://github.com/lanl-ansi/Alpine.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/lanl-ansi/Alpine.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/lanl-ansi/Alpine.jl/branch/master/badge.svg)](https://codecov.io/gh/lanl-ansi/Alpine.jl)
[![Documentation](https://github.com/lanl-ansi/Alpine.jl/actions/workflows/documentation.yml/badge.svg)](https://lanl-ansi.github.io/Alpine.jl/latest/) 
[![version](https://juliahub.com/docs/Alpine/version.svg)](https://juliahub.com/ui/Packages/Alpine/TRSJF)

ALPINE (glob(AL) o(P)timization for mixed-(I)nteger programs with (N)onlinear (E)quations), is a novel global optimization solver that uses an adaptive, piecewise convexification scheme and constraint programming methods to solve non-convex Mixed-Integer Non-Linear Programs (MINLPs) efficiently. MINLPs are typically "hard" optimization problems which appear in numerous applications (see [MINLPLib.jl](https://github.com/lanl-ansi/MINLPLib.jl)). 

Alpine is entirely built upon [JuMP](https://github.com/jump-dev/JuMP.jl) and [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) in Julia, which provides incredible flexibility for usage and further development.

Alpine globally solves a given MINLP by:

* Analyzing the problem's expressions (objective & constraints) and applies appropriate convex relaxations and polyhedral outer-approximations

* Performing sequential optimization-based bound tightening (OBBT) and an iterative MIP-based adaptive partitioning scheme via piecewise polyhedral relaxations with a guarantee of global convergence

## Installation

Install Alpine using the Julia package manager:

```julia
import Pkg
Pkg.add("Alpine")
```

## Usage with JuMP

Use Alpine with [JuMP](https://github.com/jump-dev/JuMP.jl) as follows:

```julia
using JuMP, Alpine, Ipopt, HiGHS
ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
highs = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
model = Model(
    optimizer_with_attributes(
        Alpine.Optimizer,
        "nlp_solver" => ipopt,
        "mip_solver" => highs,
    ),
)
```

## Documentation

For more details, see the [online documentation](https://lanl-ansi.github.io/Alpine.jl/latest/).

## Support problem types

Alpine can currently handle MINLPs with polynomials in constraints and/or in the objective. Currently, there is no support for exponential cones and Positive Semi-Definite (PSD) cones in MINLPs. Alpine is also a good fit for subsets of the MINLP family, for example, Mixed-Integer Quadratically Constrained Quadratic Programs (MIQCQPs), Non-Linear Programs (NLPs), etc.

For more details, check out this [video](https://www.youtube.com/watch?v=mwkhiEIS5JA) on Alpine.jl at [JuMP-dev 2018](http://www.juliaopt.org/meetings/bordeaux2018/).

## Underlying solvers

Though an MIP-based bounding algorithm implemented in Alpine is quite involved, most of the computational bottleneck arises in the underlying MIP solvers. Since every iteration of Alpine solves an MIP sub-problem, which is typically a convex MILP/MIQCQP, Alpine's run time heavily depends on the run-time of these solvers. For the best performance of Alpine, we recommend using the commercial solver [Gurobi](https://www.gurobi.com), which is available [free](https://www.gurobi.com/academia/academic-program-and-licenses/) for academic purposes. However, due to the flexibility offered by [JuMP](https://github.com/jump-dev/JuMP.jl), the following MIP and NLP solvers are supported in Alpine: 


| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)           |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)             |
| [HiGHS](https://highs.dev/)            | [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl)
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/jump-dev/Cbc.jl)                 |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/jump-dev/AmplNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/jump-dev/KNITRO.jl)           |
| [Xpress](https://www.fico.com/en/products/fico-xpress-optimization)            | [Xpress.jl](https://github.com/jump-dev/Xpress.jl)

## Bug reports and support

Please report any issues via the GitHub [issue tracker](https://github.com/lanl-ansi/Alpine.jl/issues).
All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

## Challenging Problems

We are seeking out hard benchmark instances for MINLPs. Please get in touch either by opening an issue or [privately](https://harshangrjn.github.io/#contact) if you would like to share any hard instances.

## Citing Alpine

If you find Alpine useful in your work, we kindly request that you cite the following papers ([PDF](http://harshangrjn.github.io/pdf/JOGO_2018.pdf), [PDF](http://harshangrjn.github.io/pdf/CP_2016.pdf))
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
