# Choosing Sub-Solvers

```@meta
CurrentModule = Alpine
```

The design of the AMP solver requires a variety of programming problems to be solved underneath the surface. For algorithmic performance, it's recommend that dedicated solvers to be used for these operations. The design of AMP takes advantage of MathProgBase to allow a majority of optimization softwares to be utilized easily with simple development. Currently, the following sub-solvers with Julia Interface is supported by AMP:

| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)             |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Clp.jl)                 |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)           |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/JackDunnNZ/AmplNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)           |

As the development of AMP contineous, supports fo [Mosek](http://www.mosek.com/), [GLPK](http://www.gnu.org/software/glpk/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt), [Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite) is already scheduled for the roadmap.

To use different sub-solvers, here is an example:

```julia
    using JuMP
    using Alpine
    using Gurobi, Ipopt
    m = Model()
    # Here goes the building of your model...
    setsolver(m, AlpineSolver(nlp_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(OutputFlag=0)))
```
