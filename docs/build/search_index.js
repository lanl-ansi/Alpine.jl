var documenterSearchIndex = {"docs": [

{
    "location": "intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "intro.html#POD-1",
    "page": "Introduction",
    "title": "POD",
    "category": "section",
    "text": "A Polynomial Outer-approximation Discretization Algorithm for Mixed-integer Nonlinear ProgrammingPOD.jl is a two-stage approach to strengthen piecewise convex relaxations for mixed-integer nonlinear programs (MINLP) with multi-linear terms. In the first stage, we exploit Constraint Programing techniques to contract the variable bounds. We apply feasibility-based bound contraction methods iteratively until a fixed point with respect to the bounds are achieved. In the second stage, we partition the variables domains using an adaptive multivariate partitioning scheme. Instead of equally partitioning the domains of variables appearing in multi-linear terms (predominantly common in the literature), we construct sparser partitions yet tighter relaxations by iteratively partitioning the variable domains in regions of interest (a parametrized partition around the current solution of lower-bounding MILP). This approach decouples the number of partitions from the size of the variable domains, leads to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. We further apply polyhedral cutting plane methods to handle convex relaxations of higher-order monomial terms.To use POD, external sub-solvers must be installed. See Choosing Sub-Solvers for more information."
},

{
    "location": "installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": ""
},

{
    "location": "firstdemo.html#",
    "page": "Demo",
    "title": "Demo",
    "category": "page",
    "text": ""
},

{
    "location": "firstdemo.html#Fist-Demo-1",
    "page": "Demo",
    "title": "Fist Demo",
    "category": "section",
    "text": ""
},

{
    "location": "choosingsolver.html#",
    "page": "Choosing Solvers",
    "title": "Choosing Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "choosingsolver.html#Choosing-Sub-Solvers-1",
    "page": "Choosing Solvers",
    "title": "Choosing Sub-Solvers",
    "category": "section",
    "text": "Solvers and solver-specific parameters are specified by POD, which are provided by particular solver packages. POD requires both a Non-linear Program Local Solver and a Mixed-integer Program Solver to work. For example, the Clp package exports a ClpSolver object, which can be passed to PODSolver() as follows::    using JuMP\n    using POD\n    using Gurobi, Ipopt\n    m = Model()\n    setsolver(m, PODSolver(nlp_local_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(OutputFlag=0)))"
},

]}
