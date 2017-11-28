# Choosing Sub-Solvers

The design of the POD.jl requires a variety of programming problems to be solved underneath the surface.
More importantly, the potential of POD.jl is best exploited when utilizing strong sub-solvers. The performance of the example
shown in [Usage](@ref) is approximatly 1/10 if MIP solver is Gurobi rather than Cbc.
POD.jl will customize the algorithm when utilizing different sub-solver. For example, POD.jl will utilize
tune starting values for MIP sub-solver and resolve infeasibility with local NLP/MINLP solvers.
For best algorithmic performance, it's recommend that dedicated solvers to be used for these operations.
Currently, the following sub-solvers with Julia Interface is partially supported by POD.jl:

| Solver                                                                         | Julia Package                                                |  Target Program   |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|-------------------|
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)             |   MILP, MIQCP     |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Clp.jl)                 |   MILP            |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)           |   MILP, MIQCP     |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)             |   Local NLP       |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/JackDunnNZ/AmplNLWriter.jl)   |   Local MINLP     |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)           |   Local MINLP     |

As the development of AMP contineous, supports fo [Mosek](http://www.mosek.com/), [GLPK](http://www.gnu.org/software/glpk/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), [Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite) is already scheduled for the roadmap of POD.jl.
