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
    "text": "POD.jl is a two-stage approach to strengthen piecewise convex relaxations for mixed-integer nonlinear programs (MINLP) with multi-linear terms. In the first stage, we exploit Constraint Programing techniques to contract the variable bounds. We apply feasibility-based bound contraction methods iteratively until a fixed point with respect to the bounds are achieved. In the second stage, we partition the variables domains using an adaptive multivariate partitioning scheme. Instead of equally partitioning the domains of variables appearing in multi-linear terms (predominantly common in the literature), we construct sparser partitions yet tighter relaxations by iteratively partitioning the variable domains in regions of interest (a parametrized partition around the current solution of lower-bounding MILP). This approach decouples the number of partitions from the size of the variable domains, leads to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. We further apply polyhedral cutting plane methods to handle convex relaxations of higher-order monomial terms.To use POD, external sub-solvers must be installed. See Choosing Sub-Solvers for more information."
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
    "text": "CurrentModule = PODCurrently, AMP is a private repository that is used by LANL-ANSI group, the proposed publication is unknown given the upcoming changes in JuMP and MathProgBase that AMP needs to adapt to. To install AMP,Pkg.clone(\"https://github.com/lanl-ansi/POD.git\")For developers, it is highly recommend that any further development on POD is conducted on a new branch or a forked repo."
},

{
    "location": "demo.html#",
    "page": "Demo",
    "title": "Demo",
    "category": "page",
    "text": ""
},

{
    "location": "demo.html#Demo-1",
    "page": "Demo",
    "title": "Demo",
    "category": "section",
    "text": "CurrentModule = PODConsider the following problem as an example.min_x c^Tx\nst     a_i^Tx text sense_i  b_i forall i\n         l leq x leq u"
},

{
    "location": "parameters.html#",
    "page": "Parameters",
    "title": "Parameters",
    "category": "page",
    "text": ""
},

{
    "location": "parameters.html#Parameters-1",
    "page": "Parameters",
    "title": "Parameters",
    "category": "section",
    "text": "CurrentModule = POD"
},

{
    "location": "parameters.html#General-1",
    "page": "Parameters",
    "title": "General",
    "category": "section",
    "text": "There are some general paremters that controls the behavior of the AMP:log_level: verbosity of the algorihtm (default=0), set to 1 for turning on logging, 100 for detailed debugging mode\ntimeout: total time regulated for the solver (default=Inf), in seconds\nmaxiter: total iteration allowed in the global_solve (default=999)\nrel_gap: relative gap considered for optimality during global_solve (default=1e-4), evaluated by fracUB-LBUB 100 times \ntolerance: numerical tolerance used furing the process of global_solve (default=1e-6)"
},

{
    "location": "parameters.html#Discretization-based-parameters-1",
    "page": "Parameters",
    "title": "Discretization-based parameters",
    "category": "section",
    "text": "More parameter descriptions to come..."
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
    "text": "CurrentModule = PODThe design of the AMP solver requires a variety of programming problems to be solved underneath the surface. For algorithmic performance, it's recommend that dedicated solvers to be used for these operations. The design of AMP takes advantage of MathProgBase to allow a majority of optimization softwares to be utilized easily with simple development. Currently, the following sub-solvers with Julia Interface is supported by AMP:Solver Julia Package\nCPLEX CPLEX.jl\nCbc Cbc.jl\nGurobi Gurobi.jl\nIpopt Ipopt.jl\nBonmin Bonmin.jl\nArtelys KNITRO KNITRO.jlAs the development of AMP contineous, supports fo Mosek, GLPK, NLopt, Xpress is already scheduled for the roadmap.To use different sub-solvers, here is an example:    using JuMP\n    using POD\n    using Gurobi, Ipopt\n    m = Model()\n    # Here goes the building of your model...\n    setsolver(m, PODSolver(nlp_local_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(OutputFlag=0)))"
},

{
    "location": "functions.html#",
    "page": "Methods",
    "title": "Methods",
    "category": "page",
    "text": ""
},

{
    "location": "functions.html#Functions-1",
    "page": "Methods",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = POD"
},

{
    "location": "functions.html#POD.presolve",
    "page": "Methods",
    "title": "POD.presolve",
    "category": "Function",
    "text": "presolve(m::PODNonlinearModel)\n\nAlgorithmic function that ask a loaded POD problem to perform the presolving process. It first conduct a local solve for a initial starting point followed by the bound tightening process. In cases when a local solve is infeasible, the bound tightening will continue with the bound tightening. However, as the bound tightening finishes, a second try of the local solve will be performed on a tigher domain for a initial feasible solution. The problem will be position the partitioning to the initial feasible solution or a relaxation root McCormick solution.\n\nSee Algorithm ... for details implementation.\n\n\n\n"
},

{
    "location": "functions.html#POD.global_solve",
    "page": "Methods",
    "title": "POD.global_solve",
    "category": "Function",
    "text": "global_solve(m::PODNonlinearModel)\n\nPerform the global algorithm that is based on the adaptive conexification scheme. This iterative algorithm loops over bounding_solve and local_solve for converging lower bound (relaxed problem) and upper bound (feasible problem). Each bounding_solve provides a lower bound solution that is used as a partioning point for next iteration (this feature can be modified given different add_discretization). Each local_solve provides a local serach of incumbent feasible solution. The algrithm terminates given time limits, optimality condition, or iteration limits.\n\nThe algorithm is can be reformed when add_discretization is replaced with user-defined functional input. For example, this algorithm can easily be reformed as a uniform-partitioning algorithm in other literature.\n\n\n\n"
},

{
    "location": "functions.html#POD.local_solve",
    "page": "Methods",
    "title": "POD.local_solve",
    "category": "Function",
    "text": "local_solve(m::PODNonlinearModel, presolve::Bool=false)\n\nPerform a local NLP or MINLP solve to detect feasible solution. When presolve=true, this function is utilized in the presolve for a MINLP local solve. Otherwise, a NLP problem is solved by fixing the integer variables base on the lower bound solution from bounding_solve.\n\n\n\n"
},

{
    "location": "functions.html#POD.bounding_solve",
    "page": "Methods",
    "title": "POD.bounding_solve",
    "category": "Function",
    "text": "bounding_solve(m::PODNonlinearModel; kwargs...)\n\nThis is a solving process usually deal with a MIP or MIQCP problem for lower bounds of problems. It solves the problem built upon a convexification base on a discretization Dictionary of some variables. The convexification utilized is Tighten McCormick scheme. See create_bounding_mip for more details of the problem solved here.\n\n\n\n"
},

{
    "location": "functions.html#Algorithmic-Methods-1",
    "page": "Methods",
    "title": "Algorithmic Methods",
    "category": "section",
    "text": "These are the high-level algorithmic methods:presolve\nglobal_solve\nlocal_solve\nbounding_solve"
},

{
    "location": "functions.html#Utility-Methods-1",
    "page": "Methods",
    "title": "Utility Methods",
    "category": "section",
    "text": ""
},

]}
