# Choosing Sub-Solvers

```@meta
CurrentModule = Alpine
```

Alpine's global optimization algorithm requires a few types of optimization problems to be solved efficiently. For its best performance, it is recommended that a dedicated solver is used for these purposes. The design of Alpine takes advantage of [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) to seemlessly interface with a majority of optimization software. Currently, the following sub-solvers with Julia Interface are supported by Alpine:

| Solver                                                                         | Julia Package                                                |
|--------------------------------------------------------------------------------|--------------------------------------------------------------|
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)             |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Clp.jl)                 |
| [Gurobi](http://gurobi.com/)                                                   | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)           |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)             |
| [Bonmin](https://projects.coin-or.org/Bonmin)                                  | [Bonmin.jl](https://github.com/JackDunnNZ/AmplNLWriter.jl)   |
| [Artelys KNITRO](http://artelys.com/en/optimization-tools/knitro)              | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)           |

As the development of Alpine continues, supports for solvers such as [Mosek](http://www.mosek.com/), [GLPK](http://www.gnu.org/software/glpk/), [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt) and [Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite) are in the development roadmap.

!!! tip
    Performance of Alpine is significantly faster and robust while using [Gurobi](https://www.gurobi.com) as the underlying convex mixed-integer programming (MIP) solver. Note that Gurobi's individual-usage license is available [free](https://www.gurobi.com/academia/academic-program-and-licenses/) for academic purposes. 

To solve any continuous nonlinear program to global optimality, here is an example to utilize different sub-solvers with default options set for Alpine:

```julia
using Alpine
using JuMP 
using Gurobi 
using Ipopt

# MIP optimizer
const gurobi = optimizer_with_attributes(Gurobi.Optimizer, 
                                         MOI.Silent() => true,
                                         "Presolve"   => 1) 

# NLP optimizer
const ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                        MOI.Silent() => true, 
                                        "sb" => "yes", 
                                        "max_iter"   => 9999)

# Global optimizer
const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                         "nlp_solver" => ipopt,
                                         "mip_solver" => gurobi)
m = Model(alpine)
```
Similarly, many such optimizer options and nonconvex problems (both NLPs and MINLPs) can be found in the [examples](https://github.com/lanl-ansi/Alpine.jl/tree/master/examples) folder. 