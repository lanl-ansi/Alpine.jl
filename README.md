# POD, a global MINLP solver <span style="color:black"></span>

Dev: [![Build Status](https://travis-ci.org/lanl-ansi/POD.jl.svg?branch=master)](https://travis-ci.org/lanl-ansi/POD.jl)
[![codecov](https://codecov.io/gh/lanl-ansi/POD.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/lanl-ansi/POD.jl) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://lanl-ansi.github.io/POD.jl/latest)

"Polyhedral Outer-Approximation and dynamic Discretization (POD)" is a novel global optimization algorithm that uses an adaptive convexification scheme and constraints programming methods to solve Mixed-Integer Non-Linear Programming problems (MINLPs) efficiently. MINLPs are famously known as the "hard" programming problems that exist in many applications (see this [Library](http://www.gamsworld.org/minlp/minlplib2/html/) for problem instances). POD is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadradic Convex Programming (MIQCP), Non-Linear Programming (NLP), etc.

Unlike many other state-of-the-art MINLP solvers, POD is entirely built upon [JuMP](https://github.com/JuliaOpt/JuMP.jl) and [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) Interface in Julia, which provides incredible flexibility for usage and further development.

POD solves MINLP by:

* analyzing the problem expressions (objective & constraints) to obtain a convexifyable formulation

* performing bound tightening to contract variables' feasible domains

* performing adaptive partitioning-based convexification to improve problem bounds for global convergence


## Installation

Currently, POD is a private repository that is used by the LANL-ANSI group, the proposed publication is unknown given the upcoming changes in JuMP and MathProgBase that POD needs to adapt to. To install POD,

`Pkg.clone("https://github.com/lanl-ansi/POD.git")`

For developers, it is highly recommended that any further development on POD is conducted on a new branch or a forked repo.

## Important Note
POD is a developing software and a developing research project. POD is not yet numerically robust and may require dedicated attention when encountering new instances. We embrace any issues and problems to make this open-source solver better. If you are in need of help, please contact us or post an issue.

## Citation

If you find POD useful in your work, please cite the following:
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

@techreport{bent2017polyhedral,
  title={A Polyhedral Outer-approximation, Dynamic-discretization optimization solver, 1. x},
  author={Bent, Rusell and Nagarajan, Harsha and Sundar, Kaarthik and Wang, Site and Hijazi, Hassan},
  year={2017},
  institution={Los Alamos National Laboratory (LANL), Los Alamos, NM (United States)}
}
```
