var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "#Alpine-1",
    "page": "Introduction",
    "title": "Alpine",
    "category": "section",
    "text": "Alpine.jl is a two-stage approach to strengthen piecewise convex relaxations for mixed-integer nonlinear programs (MINLP) with multi-linear terms. In the first stage, we exploit Constraint Programing techniques to contract the variable bounds. We apply feasibility-based bound contraction methods iteratively until a fixed point with respect to the bounds are achieved. In the second stage, we partition the variables domains using an adaptive multivariate partitioning scheme. Instead of equally partitioning the domains of variables appearing in multi-linear terms (predominantly common in the literature), we construct sparser partitions yet tighter relaxations by iteratively partitioning the variable domains in regions of interest (a parametrized partition around the current solution of lower-bounding MILP). This approach decouples the number of partitions from the size of the variable domains, leads to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. We further apply polyhedral cutting plane methods to handle convex relaxations of higher-order monomial terms.To use Alpine, external sub-solvers must be installed. See Choosing Sub-Solvers for more information."
},

{
    "location": "installation/#",
    "page": "How to Use",
    "title": "How to Use",
    "category": "page",
    "text": ""
},

{
    "location": "installation/#Installation-1",
    "page": "How to Use",
    "title": "Installation",
    "category": "section",
    "text": "CurrentModule = AlpineCurrently, AMP is a private repository that is used by LANL-ANSI group, the proposed publication is unknown given the upcoming changes in JuMP and MathProgBase that AMP needs to adapt to. To install AMP,Pkg.clone(\"https://github.com/lanl-ansi/Alpine.git\")For developers, it is highly recommend that any further development on Alpine is conducted on a new branch or a forked repo."
},

{
    "location": "choosingsolver/#",
    "page": "Choosing Sub-Solvers",
    "title": "Choosing Sub-Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "choosingsolver/#Choosing-Sub-Solvers-1",
    "page": "Choosing Sub-Solvers",
    "title": "Choosing Sub-Solvers",
    "category": "section",
    "text": "CurrentModule = AlpineThe design of the AMP solver requires a variety of programming problems to be solved underneath the surface. For algorithmic performance, it\'s recommend that dedicated solvers to be used for these operations. The design of AMP takes advantage of MathProgBase to allow a majority of optimization softwares to be utilized easily with simple development. Currently, the following sub-solvers with Julia Interface is supported by AMP:Solver Julia Package\nCPLEX CPLEX.jl\nCbc Cbc.jl\nGurobi Gurobi.jl\nIpopt Ipopt.jl\nBonmin Bonmin.jl\nArtelys KNITRO KNITRO.jlAs the development of AMP contineous, supports fo Mosek, GLPK, NLopt, Xpress is already scheduled for the roadmap.To use different sub-solvers, here is an example:    using JuMP\n    using Alpine\n    using Gurobi, Ipopt\n    m = Model()\n    # Here goes the building of your model...\n    setsolver(m, AlpineSolver(nlp_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(OutputFlag=0)))"
},

{
    "location": "algorithm/#",
    "page": "Algorithm",
    "title": "Algorithm",
    "category": "page",
    "text": ""
},

{
    "location": "algorithm/#Demo-1",
    "page": "Algorithm",
    "title": "Demo",
    "category": "section",
    "text": "CurrentModule = AlpineConsider the following problem as an example.min_x c^Tx\nst     a_i^Tx text sense_i  b_i forall i\n         l leq x leq u"
},

{
    "location": "expression/#",
    "page": "Expression Guideline",
    "title": "Expression Guideline",
    "category": "page",
    "text": ""
},

{
    "location": "parameters/#",
    "page": "Parameters",
    "title": "Parameters",
    "category": "page",
    "text": ""
},

{
    "location": "parameters/#Parameters-1",
    "page": "Parameters",
    "title": "Parameters",
    "category": "section",
    "text": "CurrentModule = Alpine"
},

{
    "location": "parameters/#General-1",
    "page": "Parameters",
    "title": "General",
    "category": "section",
    "text": "There are some general parameters that controls the behavior of the AMP:loglevel(default=0): verbosity of the algorithm, set to 1 for turning on logging, 100 for detailed debugging mode\ntimeout(default=Inf): total time regulated for the solver in seconds\nmaxiter(default=999): total iteration allowed in the global_solve\nrelgap(default=1e-4): relative gap considered for optimality during global_solve. Bounds are evaluated using fracUB-LBUB 100 times \ntol(default=1e-6): numerical tol used during the process of global_solve"
},

{
    "location": "parameters/#Discretization-based-parameters-1",
    "page": "Parameters",
    "title": "Discretization-based parameters",
    "category": "section",
    "text": "disc_ratio(default=4): used during add_adpative_partition for measuring the radius of new partitioning relative to the active domain\ndisc_var_pick(default=0): controls algorithm/methods used for selecting variables for discretization, 0 is for max-cover, 1 is for minimum-vertex-cover. This parameter allows functional inputs.\ndisc_add_partition_method: allows functional input on how new partitions should be constructed. This parameter is not fully stable."
},

{
    "location": "parameters/#Presolve-Parameters-1",
    "page": "Parameters",
    "title": "Presolve Parameters",
    "category": "section",
    "text": "presolve_track_time(default=false): consider presolve time as the total time or not\npresolve_bt(default=false): perform built-in bound tightening presolve procedure\npresolve_maxiter(default=9999): maximum iteration allowed using presolve process\npresolve_bt_width_tol(default=1e-3): independent numerical tol used in presolve for bound tightening procedure. Note that this procedure is more sensitive to the tol in here. Small tol is more likely to results in strange presolve behavior.\npresolve_bt_algo(default=1): method used to do built-in bound tightening, choose 1 for regular bounding tightening,  2 for Tighten McCormick bound tightening.\npresolve_bt_relax(default=false): whether to relax the bounding tightening MILP solved or not\npresolve_bt_mip_timeout(default=Inf): time limit used for individual MILP solved during bound tightening presolveMore parameter descriptions to come..."
},

{
    "location": "functions/#",
    "page": "Methods",
    "title": "Methods",
    "category": "page",
    "text": ""
},

{
    "location": "functions/#Functions-1",
    "page": "Methods",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = Alpine"
},

{
    "location": "functions/#Alpine.presolve",
    "page": "Methods",
    "title": "Alpine.presolve",
    "category": "function",
    "text": "presolve(m::AlpineNonlinearModel)\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.global_solve",
    "page": "Methods",
    "title": "Alpine.global_solve",
    "category": "function",
    "text": "global_solve(m::AlpineNonlinearModel)\n\nPerform global optimization algorithm that is based on the adaptive piecewise convexification. This iterative algorithm loops over bounding_solve and local_solve until the optimality gap between the lower bound (relaxed problem with min. objective) and the upper bound (feasible problem) is within the user prescribed limits. Each bounding_solve provides a lower bound that serves as the partioning point for the next iteration (this feature can be modified given a different add_adaptive_partition). Each local_solve provides an incumbent feasible solution. The algorithm terminates when atleast one of these conditions are satisfied: time limit, optimality condition, or iteration limit.\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.local_solve",
    "page": "Methods",
    "title": "Alpine.local_solve",
    "category": "function",
    "text": "local_solve(m::AlpineNonlinearModel, presolve::Bool=false)\n\nPerform a local NLP or MINLP solve to obtain a feasible solution. The presolve option is set to true when the function is invoked in presolve. Otherwise, the function is invoked from bounding_solve.\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.bounding_solve",
    "page": "Methods",
    "title": "Alpine.bounding_solve",
    "category": "function",
    "text": "bounding_solve(m::AlpineNonlinearModel; kwargs...)\n\nThis is a solving process usually deal with a MIP or MIQCP problem for lower bounds of problems. It solves the problem built upon a convexification base on a discretization Dictionary of some variables. The convexification utilized is Tighten McCormick scheme. See create_bounding_mip for more details of the problem solved here.\n\n\n\n\n\n"
},

{
    "location": "functions/#High-level-Algorithmic-Operations-1",
    "page": "Methods",
    "title": "High-level Algorithmic Operations",
    "category": "section",
    "text": "These are the high-level algorithmic methods:presolve\nglobal_solve\nlocal_solve\nbounding_solve"
},

{
    "location": "functions/#Adapative-Partitioning-Methods-1",
    "page": "Methods",
    "title": "Adapative Partitioning Methods",
    "category": "section",
    "text": "create_bounding_mip\npick_disc_vars\nfix_domains\nmin_vertex_cover\nmax_cover"
},

{
    "location": "functions/#Alpine.bound_tightening",
    "page": "Methods",
    "title": "Alpine.bound_tightening",
    "category": "function",
    "text": "bound_tightening(m::AlpineNonlinearModel)\n\nEntry point for the bound-tightening algorithm. The aim of the bound-tightening algorithm is to tighten the variable bounds, if possible.\n\nCurrently, two bounding tightening method is implemented minmax_bound_tightening.\n\n* Bound-tightening with basic McCormick\n* Bound-tightening with McCormick partitions: (3 partitions around the local feasible solution)\nIf no local feasible solution is obtained, the algorithm defaults to bound-tightening with basic McCormick\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.minmax_bound_tightening",
    "page": "Methods",
    "title": "Alpine.minmax_bound_tightening",
    "category": "function",
    "text": "minmax_bound_tightening(m:AlpineNonlinearModel; use_bound::Bool=true, use_tmc::Bool)\n\nThis function implements the bound-tightening algorithm to tighten the variable bounds. It utilizes either the basic McCormick relaxation or the Tightened McCormick relaxation (TMC) to tighten the bounds. The TMC has additional binary variables for partitioning.\n\nThe algorithm as two main parameters. The first is the use_tmc, which when set to true invokes the algorithm on the TMC relaxation. The second parameter use_bound takes in the objective value of the local solve solution stored in best_sol. The use_bound option is set to true when the local solve is successful is obtaining a feasible solution and the bound-tightening is to be performed using the objective value of the feasible solution, else this parameter is set to false\n\nSeveral other parameters are available for the presolve algorithm tuning. For more details, see Parameters.\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.create_bound_tightening_model",
    "page": "Methods",
    "title": "Alpine.create_bound_tightening_model",
    "category": "function",
    "text": "create_bound_tightening_model(m::AlpineNonlinearModel, discretization::Dict, bound::Float64)\n\nThis function takes in the initial discretization information and builds a bound-tightening model. It is an algorithm specific function called by minmax_bound_tightening\n\n\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.solve_bound_tightening_model",
    "page": "Methods",
    "title": "Alpine.solve_bound_tightening_model",
    "category": "function",
    "text": "solve_bound_tightening_model(m::AlpineNonlinearModel)\n\nA function that solves the min and max bound-tightening model.\n\n\n\n\n\n"
},

{
    "location": "functions/#Alpine.resolve_var_bounds",
    "page": "Methods",
    "title": "Alpine.resolve_var_bounds",
    "category": "function",
    "text": "resolve_var_bounds(m::AlpineNonlinearModel)\n\nResolve the bounds of the lifted variable using the information in lvartight and uvartight. This method only takes in known or trivial bounds information to reason lifted variable bound to avoid the cases of infinity bounds.\n\n\n\n\n\nresolve_var_bounds(nonconvex_terms::Dict, discretization::Dict)\n\nFor discretization to be performed, we do not allow for a variable being discretized to have infinite bounds.\nThe lifted variables will have infinite bounds and the function infers bounds on these variables. This process\ncan help speed up the subsequent solve in subsequent iterations.\n\nOnly used in presolve bound tightening\n\n\n\n\n\n"
},

{
    "location": "functions/#Presolve-Methods-1",
    "page": "Methods",
    "title": "Presolve Methods",
    "category": "section",
    "text": "bound_tightening\nminmax_bound_tightening\ncreate_bound_tightening_model\nsolve_bound_tightening_model\nresolve_var_bounds"
},

{
    "location": "functions/#Utility-Methods-1",
    "page": "Methods",
    "title": "Utility Methods",
    "category": "section",
    "text": "update_var_bounds\ndiscretization_to_bounds\ninit_disc\nto_discretization\nflatten_discretization\nadd_adpative_partition\nupdate_mip_time_limit\nfetch_timeleft_symbol"
},

{
    "location": "hacking/#",
    "page": "Hacking Alpine",
    "title": "Hacking Alpine",
    "category": "page",
    "text": ""
},

{
    "location": "hacking/#Hacking-Solver-1",
    "page": "Hacking Alpine",
    "title": "Hacking-Solver",
    "category": "section",
    "text": "CurrentModule = AlpineThis page give more detailed tutorial on how to hack this solver for different behaviors..."
},

]}
