# POD, A Global Solver for Nonconvex MINLPs <span style="color:black"></span>

Dev: [![Build Status](https://travis-ci.org/lanl-ansi/POD.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/POD.jl)
[![codecov](https://codecov.io/gh/lanl-ansi/POD.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/POD.jl)

"POD: Piecewise convex relaxation (P), Outer-approximation (O), and Dynamic discretization (D)" is a novel global optimization algorithm that uses an adaptive convexification scheme and constraints programming methods to solve Mixed-Integer Non-Linear Programs (non-convex MINLPs) efficiently. MINLPs are famously known as the "hard" programming problems that exist in many applications (see this [MINLPLibJuMP.jl](https://github.com/lanl-ansi/MINLPLibJuMP.jl) for problem instances). POD is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadradic Convex Programming (MIQCP), Non-Linear Programming (NLP), etc.

Unlike many other state-of-the-art MINLP solvers, POD is entirely built upon [JuMP](https://github.com/JuliaOpt/JuMP.jl) and [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) Interface in Julia, which provides incredible flexibility for usage and further development.

POD globally solves a given MINLP by:

* Analyzing the problem's expressions (objective & constraints) and applies approporite convex relaxations

* Performing novel adaptive/dynamic partitioning methods to create piecewise relaxations, bound tightening and polyhedral outer-approximations to guarantee global convergence

 **Illustration of POD's dynamic partitioning and outer-approximation on simple functions** ([Source](https://arxiv.org/abs/1707.02514))
 
<p align="center"> <img src="https://github.com/lanl-ansi/POD.jl/blob/master/Dynamic_partitions_github.png" width="580" class="centerImage"> </p>

{{< youtube mwkhiEIS5JA >}}
	 


## Installation

POD is a currently a unregistered package, with it's repository under the LANL-ANSI group, which can be installed through the Julia package manager:

`Pkg.clone("https://github.com/lanl-ansi/POD.git")`

Developers: Any further development of POD can be conducted on a new branch or a forked repo.

## Underlying solvers

Though the algorithm implemented in POD is quite involved, most of the hard work and computational bottleneck would arise in the underlying solvers. Since every iteration of POD solves a subproblem to optimality, which is typically a convex MILP/MIQCQP and solves a nonconvex NLP/MINLP to local optimality, POD's run time heavily depends on the quality of these solvers. For best performance of POD, use commercial solvers such as CPLEX/Gurobi. However, due to the flexibility offered by JuMP/MathProgBase, the following solvers are supported in POD: 


| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)             |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Clp.jl)                 |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)           |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/JackDunnNZ/PODlNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)           |

As the development of POD continues, support for solvers like [Mosek](http://www.mosek.com/), [Pajarito](https://github.com/JuliaOpt/Pajarito.jl), [GLPK](http://www.gnu.org/software/glpk/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), [Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite) is in the roadmap.

## Bug reports and support
Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

[issue tracker]: https://github.com/lanl-ansi/POD.jl/issues


## Citing POD

If you find POD useful in your work, we kindly request that you cite the following papers
```
@article{nagarajan2017adaptive,
  title={An Adaptive, Multivariate Partitioning Algorithm for Global Optimization of Nonconvex Programs},
  author={Nagarajan, Harsha and Lu, Mowen and Wang, Site and Bent, Russell and Sundar, Kaarthik},
  journal={arXiv preprint arXiv:1707.02514},
  year={2017}
}

@inproceedings{nagarajan2016tightening,
  title={Tightening {McC}ormick relaxations for nonlinear programs via dynamic multivariate partitioning},
  author={Nagarajan, Harsha and Lu, Mowen and Yamangil, Emre and Bent, Russell},
  booktitle={International Conference on Principles and Practice of Constraint Programming},
  pages={369--387},
  year={2016},
  organization={Springer}
}
```
