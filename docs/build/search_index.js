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
    "text": "CurrentModule = PODConsider the following problem as an example.min_x c^Tx\nst     a_i^Tx text sense_i  b_i forall i\n         l leq x leq u"
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
    "text": "There are some general paremters that controls the behavior of the AMP:log_level(default=0): verbosity of the algorihtm, set to 1 for turning on logging, 100 for detailed debugging mode\ntimeout(default=Inf): total time regulated for the solver in seconds\nmaxiter(default=999): total iteration allowed in the global_solve\nrel_gap(default=1e-4): relative gap considered for optimality during global_solve. Bounds are evaluated using fracUB-LBUB 100 times \ntol(default=1e-6): numerical tol used furing the process of global_solve"
},

{
    "location": "parameters.html#Discretization-based-parameters-1",
    "page": "Parameters",
    "title": "Discretization-based parameters",
    "category": "section",
    "text": "discretization_ratio(default=4): used during add_discretization for measuring the radius of new partitioning relative to the active domain\ndiscretization_var_pick_algo(default=0): controls algorithm/methods used for selecting variables for discretization, 0 is for max-cover, 1 is for minimum-vertex-cover. This parameter allows functional inputs.\ndiscretization_add_partition_method: allows functional input on how new partitions should be constructed. This paremeter is not fully stable."
},

{
    "location": "parameters.html#Presolve-Parameters-1",
    "page": "Parameters",
    "title": "Presolve Parameters",
    "category": "section",
    "text": "presolve_track_time(default=false): consier presolve time as the total time or not\npresolve_perform_bound_tightening(default=false): perform built-in bound tightening presolve procedure\npresolve_maxiter(default=9999): maximum iteration allowed using presolve process\npresolve_bt_width_tol(default=1e-3): independent numerical tol used in presolve for bound tightening procedure. Note that this procedure is more sensitive to the tol in here. Small tol is more likely to results in strange presolve behvaior.\npresolve_bound_tightening_algo(default=1): method used to do built-in bound tightening, choose 1 for regular bounding tightening,  2 for Tighten McCormick bound tightening.\npresolve_mip_relaxation(default=false): whether to relax the bounding tightening MILP solved or not\npresolve_mip_timelimit(default=Inf): time limit used for invidiual MILP solved during bound tightening presolveMore parameter descriptions to come..."
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
    "text": "CurrentModule = PODThe design of the AMP solver requires a variety of programming problems to be solved underneath the surface. For algorithmic performance, it's recommend that dedicated solvers to be used for these operations. The design of AMP takes advantage of MathProgBase to allow a majority of optimization softwares to be utilized easily with simple development. Currently, the following sub-solvers with Julia Interface is supported by AMP:Solver Julia Package\nCPLEX CPLEX.jl\nCbc Cbc.jl\nGurobi Gurobi.jl\nIpopt Ipopt.jl\nBonmin Bonmin.jl\nArtelys KNITRO KNITRO.jlAs the development of AMP contineous, supports fo Mosek, GLPK, NLopt, Xpress is already scheduled for the roadmap.To use different sub-solvers, here is an example:    using JuMP\n    using POD\n    using Gurobi, Ipopt\n    m = Model()\n    # Here goes the building of your model...\n    setsolver(m, PODSolver(nlp_local_solver=IpoptSolver(print_level=0), mip_solver=GurobiSolver(OutputFlag=0)))"
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
    "text": "presolve(m::PODNonlinearModel)\n\nFunction that perfoms a presolve on the user-supplied nonlinear program. The presolve first provides the model to a local solver to obtain an initial feasible solution. If the local solver returns a feasible solution then the objective value of the feasible solution is used in conjunction with the original model for a bound-tightening procedure. If the local solver reports infeasibililty the bound-tightening procedure does not use the objective value. Furthermore, in this case, a second local solve is attempted using the tightened bounds.\n\nIf a local solution is not obtained eved after the second solve then an initial McCormick solve is performed. The local solution (if available) or the initial McCormick solution (if infeasible after two local solve tries) is then used to partition the variables for the subsequent Adaptive Multivariate Partitioning algorithm iterations.\n\n\n\n"
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
    "text": "local_solve(m::PODNonlinearModel, presolve::Bool=false)\n\nPerform a local NLP or MINLP solve to obtain a feasible solution. The presolve option is set to true when the function is invoked in presolve. Otherwise, the function is invoked from bounding_solve.\n\n\n\n"
},

{
    "location": "functions.html#POD.bounding_solve",
    "page": "Methods",
    "title": "POD.bounding_solve",
    "category": "Function",
    "text": "bounding_solve(m::PODNonlinearModel; kwargs...)\n\nThis is a solving process usually deal with a MIP or MIQCP problem for lower bounds of problems. It solves the problem built upon a convexification base on a discretization Dictionary of some variables. The convexification utilized is Tighten McCormick scheme. See create_bounding_mip for more details of the problem solved here.\n\n\n\n"
},

{
    "location": "functions.html#High-level-Algorithmic-Operations-1",
    "page": "Methods",
    "title": "High-level Algorithmic Operations",
    "category": "section",
    "text": "These are the high-level algorithmic methods:presolve\nglobal_solve\nlocal_solve\nbounding_solve"
},

{
    "location": "functions.html#POD.create_bounding_mip",
    "page": "Methods",
    "title": "POD.create_bounding_mip",
    "category": "Function",
    "text": "create_bounding_mip(m::PODNonlinearModel; use_discretization::Dict)\n\nSet up a JuMP MILP bounding model base on variable domain partitioning information stored in use_discretization. By default, if use_discretization is not provided, it will use m.discretizations store in the POD model. The basic idea of this MILP bounding model is to use Tighten McCormick to convexify the original Non-convex region. Among all presented partitionings, the bounding model will choose one specific partition as the lower bound solution. The more partitions there are, the better or finer bounding model relax the original MINLP while the more efforts required to solve this MILP is required.\n\nThis function is implemented in the following manner:\n\n* [`post_amp_vars`](@ref): post original and lifted variables\n* [`post_amp_lifted_constraints`](@ref): post original and lifted constraints\n* [`post_amp_lifted_obj`](@ref): post original or lifted objective function\n* [`post_amp_mccormick`](@ref): post Tighen McCormick variables and constraints base on `discretization` information\n\nMore specifically, the Tightening McCormick used here can be genealized in the following mathematcial formulation. Consider a nonlinear term\n\nbeginsubequations\nbeginalign\n   widehatx_ij geq (mathbfx_i^lcdothatmathbfy_i) x_j + (mathbfx_j^lcdothatmathbfy_j) x_i - (mathbfx_i^lcdothatmathbfy_i)(mathbfx_j^lcdothatmathbfy_j) \n   widehatx_ij geq (mathbfx_i^ucdothatmathbfy_i) x_j + (mathbfx_j^ucdothatmathbfy_j) x_i - (mathbfx_i^ucdothatmathbfy_i)(mathbfx_j^ucdothatmathbfy_j) \n   widehatx_ij leq (mathbfx_i^lcdothatmathbfy_i) x_j + (mathbfx_j^ucdothatmathbfy_j) x_i - (mathbfx_i^lcdothatmathbfy_i)(mathbfx_j^ucdothatmathbfy_j) \n   widehatx_ij leq (mathbfx_i^ucdothatmathbfy_i) x_j + (mathbfx_j^lcdothatmathbfy_j) x_i - (mathbfx_i^ucdothatmathbfy_i)(mathbfx_j^lcdothatmathbfy_j) \n    mathbfx_i^ucdothatmathbfy_i) geq x_i geq mathbfx_i^lcdothatmathbfy_i) \n    mathbfx_j^ucdothatmathbfy_j) geq x_j geq mathbfx_j^lcdothatmathbfy_j) \n   sum hatmathbfy_i = 1   sum hatmathbfy_j_k = 1 \n   hatmathbfy_i in 01 hatmathbfy_j in 01\nendalign\nendsubequations\n\n\n\n"
},

{
    "location": "functions.html#POD.pick_vars_discretization",
    "page": "Methods",
    "title": "POD.pick_vars_discretization",
    "category": "Function",
    "text": "pick_vars_discretization(m::PODNonlinearModel)\n\nThis function helps pick the variables for discretization. The method chosen depends on user-inputs. In case when indices::Int is provided, the method is chosen as built-in method. Currently, there exist two built-in method:\n\n* `max-cover(m.discretization_var_pick_algo=0, default)`: pick all variables involved in the non-linear term for discretization\n* `min-vertex-cover(m.discretization_var_pick_algo=1)`: pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified\n\nFor advance usage, m.discretization_var_pick_algo allows ::Function inputs. User is required to perform flexible methods in choosing the non-linear variable. For more information, read more details at Hacking Solver.\n\n\n\n"
},

{
    "location": "functions.html#POD.fix_domains",
    "page": "Methods",
    "title": "POD.fix_domains",
    "category": "Function",
    "text": "fix_domains(m::PODNonlinearModel)\n\nThis function is used to fix variables to certain domains during the local solve process in the global_solve. More specifically, it is used in local_solve to fix binary and integer variables to lower bound solutions and discretizing varibles to the active domain according to lower bound solution.\n\n\n\n"
},

{
    "location": "functions.html#POD.min_vertex_cover",
    "page": "Methods",
    "title": "POD.min_vertex_cover",
    "category": "Function",
    "text": "min_vertex_cover(m:PODNonlinearModel)\n\nA built-in method for selecting variables for discretization.\n\n\n\n"
},

{
    "location": "functions.html#POD.max_cover",
    "page": "Methods",
    "title": "POD.max_cover",
    "category": "Function",
    "text": "max_cover(m:PODNonlinearModel)\n\nA built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.\n\n\n\n"
},

{
    "location": "functions.html#Adapative-Partitioning-Methods-1",
    "page": "Methods",
    "title": "Adapative Partitioning Methods",
    "category": "section",
    "text": "create_bounding_mip\npick_vars_discretization\nfix_domains\nmin_vertex_cover\nmax_cover"
},

{
    "location": "functions.html#POD.bound_tightening",
    "page": "Methods",
    "title": "POD.bound_tightening",
    "category": "Function",
    "text": "bound_tightening(m::PODNonlinearModel)\n\nEntry point for the bound-tightening algorithm. The aim of the bound-tightening algorithm is to tighten the variable bounds, if possible.\n\nCurrently, two bounding tightening method is implemented minmax_bound_tightening.\n\n* Bound-tightening with basic McCormick\n* Bound-tightening with McCormick partitions: (3 partitions around the local feasible solution)\nIf no local feasible solution is obtained, the algorithm defaults to bound-tightening with basic McCormick\n\n\n\n"
},

{
    "location": "functions.html#POD.minmax_bound_tightening",
    "page": "Methods",
    "title": "POD.minmax_bound_tightening",
    "category": "Function",
    "text": "minmax_bound_tightening(m:PODNonlinearModel; use_bound::Bool=true, use_tmc::Bool)\n\nThis function implements the bound-tightening algorithm to tighten the variable bounds. It utilizes either the basic McCormick relaxation or the Tightened McCormick relaxation (TMC) to tighten the bounds. The TMC has additional binary variables for partitioning.\n\nThe algorithm as two main parameters. The first is the use_tmc, which when set to true invokes the algorithm on the TMC relaxation. The second parameter use_bound takes in the objective value of the local solve solution stored in best_sol. The use_bound option is set to true when the local solve is successful is obtaining a feasible solution and the bound-tightening is to be performed using the objective value of the feasible solution, else this parameter is set to false\n\nSeveral other parameters are available for the presolve algorithm tuning. For more details, see Parameters.\n\n\n\n"
},

{
    "location": "functions.html#POD.create_bound_tightening_model",
    "page": "Methods",
    "title": "POD.create_bound_tightening_model",
    "category": "Function",
    "text": "create_bound_tightening_model(m::PODNonlinearModel, discretization::Dict, bound::Float64)\n\nThis function takes in the initial discretization information and builds a bound-tightening model. It is an algorithm specific function called by minmax_bound_tightening\n\n\n\n\n\n"
},

{
    "location": "functions.html#POD.solve_bound_tightening_model",
    "page": "Methods",
    "title": "POD.solve_bound_tightening_model",
    "category": "Function",
    "text": "solve_bound_tightening_model(m::PODNonlinearModels)\n\nA function that solves the min and max bound-tightening model.\n\n\n\n"
},

{
    "location": "functions.html#POD.resolve_lifted_var_bounds",
    "page": "Methods",
    "title": "POD.resolve_lifted_var_bounds",
    "category": "Function",
    "text": "resolve_lifted_var_bounds(nonlinear_info::Dict, discretization::Dict)\n\nFor discretization to be performed, we do not allow for a variable being discretized to have infinite bounds. The lifted variables will have infinite bounds and the function infers bounds on these variables. This process can help speed up the subsequent solve in subsequent iterations.\n\n\n\n"
},

{
    "location": "functions.html#Presolve-Methods-1",
    "page": "Methods",
    "title": "Presolve Methods",
    "category": "section",
    "text": "bound_tightening\nminmax_bound_tightening\ncreate_bound_tightening_model\nsolve_bound_tightening_model\nresolve_lifted_var_bounds"
},

{
    "location": "functions.html#POD.update_var_bounds",
    "page": "Methods",
    "title": "POD.update_var_bounds",
    "category": "Function",
    "text": "update_var_bounds(m::PODNonlinearModel, discretization::Dict; len::Float64=length(keys(discretization)))\n\nThis function take in a dictionary-based discretization information and convert them into two bounds vectors (l_var, u_var) by picking the smallest and largest numbers. User can specify a certain length that may contains variables that is out of the scope of discretization.\n\nOutput::\n\nl_var::Vector{Float64}, u_var::Vector{Float64}\n\n\n\n"
},

{
    "location": "functions.html#POD.discretization_to_bounds",
    "page": "Methods",
    "title": "POD.discretization_to_bounds",
    "category": "Function",
    "text": "discretization_to_bounds(d::Dict, l::Int)\n\nSame as update_var_bounds\n\n\n\n"
},

{
    "location": "functions.html#POD.initialize_discretization",
    "page": "Methods",
    "title": "POD.initialize_discretization",
    "category": "Function",
    "text": "initialize_discretization(m::PODNonlinearModel)\n\nThis function initialize the dynamic discretization used for any bounding models. By default, it takes (.l_var_orig, .u_var_orig) as the base information. User is allowed to use alternative bounds for initializing the discretization dictionary. The output is a dictionary with MathProgBase variable indices keys attached to the :PODNonlinearModel.discretization.\n\n\n\n"
},

{
    "location": "functions.html#POD.to_discretization",
    "page": "Methods",
    "title": "POD.to_discretization",
    "category": "Function",
    "text": "to_discretization(m::PODNonlinearModel, lbs::Vector{Float64}, ubs::Vector{Float64})\n\nUtility functions to convert bounds vectors to Dictionary based structures that is more suitable for partition operations.\n\n\n\n"
},

{
    "location": "functions.html#POD.flatten_discretization",
    "page": "Methods",
    "title": "POD.flatten_discretization",
    "category": "Function",
    "text": "flatten_discretization(discretization::Dict)\n\nUtility functions to eliminate all partition on discretizing variable and keep the loose bounds.\n\n\n\n"
},

{
    "location": "functions.html#POD.add_discretization",
    "page": "Methods",
    "title": "POD.add_discretization",
    "category": "Function",
    "text": "add_discretization(m::PODNonlinearModel; use_discretization::Dict, use_solution::Vector)\n\nBasic built-in method used to add a new partition on feasible domains of discretizing variables. This method make modification in .discretization\n\nConsider original partition [0, 3, 7, 9], where LB/any solution is 4. Use ^ as the new partition, \"|\" as the original partition\n\nA case when discretize ratio = 4 | –––– | - ^ – * – ^ –– | –––– | 0          3  3.5   4   4.5     7          9\n\nA special case when discretize ratio = 2 | –––– | –– * –– ^ –– | –––– | 0          3      4      5      7          9\n\nThere are two options for this function,\n\n* `use_discretization(default=m.discretization)`:: to regulate which is the base to add new partitions on\n* `use_solution(default=m.best_bound_sol)`:: to regulate which solution to use when adding new partitions on\n\nThis function belongs to the hackable group, which means it can be replaced by the user to change the behvaior of the solver.\n\n\n\n"
},

{
    "location": "functions.html#POD.update_mip_time_limit",
    "page": "Methods",
    "title": "POD.update_mip_time_limit",
    "category": "Function",
    "text": "update_mip_time_limit(m::PODNonlinearModel)\n\nAn utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.\n\n\n\n"
},

{
    "location": "functions.html#POD.fetch_timeleft_symbol",
    "page": "Methods",
    "title": "POD.fetch_timeleft_symbol",
    "category": "Function",
    "text": "fetch_timeleft_symbol(m::PODNonlinearModel)\n\nAn utility function used to recognize differnt sub-solvers return the timelimit setup keywords.\n\n\n\n"
},

{
    "location": "functions.html#Utility-Methods-1",
    "page": "Methods",
    "title": "Utility Methods",
    "category": "section",
    "text": "update_var_bounds\ndiscretization_to_bounds\ninitialize_discretization\nto_discretization\nflatten_discretization\nadd_discretization\nupdate_mip_time_limit\nfetch_timeleft_symbol"
},

{
    "location": "hacking.html#",
    "page": "Hacking",
    "title": "Hacking",
    "category": "page",
    "text": ""
},

{
    "location": "hacking.html#Hacking-Solver-1",
    "page": "Hacking",
    "title": "Hacking-Solver",
    "category": "section",
    "text": "CurrentModule = PODThis page give more detailed tutorial on how to hack this solver for different behaviors..."
},

]}
