# AMP(POD), a novel MINLP solver

Adaptive Multi-variant Partitioning (AMP/POD) is a novel global algorithm that use adaptive convexification scheme and constraints programming methods to solve mixed-integer non-linear programming problems (MINLPs) efficiently. MINLPs are famously known as the "hard" programming problems that exist in many applications (see this [Library](http://www.gamsworld.org/minlp/minlplib2/html/) for problem instances). AMP is also a fit for many relaxations of the MINLPs, such as Mixed-Integer Quadradic Convex Programming (MIQCP), Non-Linear Programming (NLP), etc.

Unlike many other state-of-the-art MINLP solvers, AMP is entirely built upon [JuMP](https://github.com/JuliaOpt/JuMP.jl) and [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) Interface in Julia, which provides incredible flexibility for usage and further development.

AMP solves MINLP by:

* analyzing the problem expressions (objective & constraints) to obtain a convexify-able formulation

* perform bound tightening algorithm to contract variable feasible domains
    
* perform adaptive partitioning-based convexification to improve problem bounds for global convergence


# Installation

Currently, AMP is a private repository that is used by LANL-ANSI group, the proposed publication is unknown given the upcoming changes in JuMP and MathProgBase that AMP needs to adapt to. To install AMP,

`Pkg.clone("https://github.com/lanl-ansi/POD.git")`

For developers, it is highly recommend that any further development on POD is conducted on a new branch or a forked repo.

# Usage

The design of the AMP solver requires a variety of programming problems to be solved underneath the surface. For algorithmic performance, it's recommend that dedicated solvers to be used for these operations. The design of AMP takes advantage of MathProgBase to allow a majority of optimization softwares to be utilized easily with simple development. Currently, the following sub-solvers with Julia Interface is supported by AMP:

| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)             |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Clp.jl)                 |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)           |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/JackDunnNZ/AmplNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)           |

As the development of AMP contineous, supports fo [Mosek](http://www.mosek.com/), [Pajarito](https://github.com/JuliaOpt/Pajarito.jl), [GLPK](http://www.gnu.org/software/glpk/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), [Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite) is already scheduled for the roadmap.

# Citation

If you find AMP useful in your work, we kindly request your citation.
