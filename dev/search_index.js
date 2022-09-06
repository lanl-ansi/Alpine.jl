var documenterSearchIndex = {"docs":
[{"location":"choosingsolver/#Choosing-Sub-Solvers","page":"Choosing Sub-solvers","title":"Choosing Sub-Solvers","text":"","category":"section"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"CurrentModule = Alpine","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"Alpine's global optimization algorithm requires a few types of optimization problems to be solved efficiently. For its best performance, it is recommended that a dedicated solver is used for these purposes. The design of Alpine takes advantage of MathOptInterface to seemlessly interface with a majority of optimization software. Currently, the following sub-solvers with Julia Interface are supported by Alpine:","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"Solver Julia Package\nCPLEX CPLEX.jl\nCbc Cbc.jl\nGurobi Gurobi.jl\nIpopt Ipopt.jl\nBonmin Bonmin.jl\nArtelys KNITRO KNITRO.jl","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"As the development of Alpine continues, supports for solvers such as Mosek, GLPK, NLopt and Xpress are in the development roadmap.","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"tip: Tip\nPerformance of Alpine is significantly faster and robust while using Gurobi as the underlying convex mixed-integer programming (MIP) solver. Note that Gurobi's individual-usage license is available free for academic purposes. ","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"To solve any continuous nonlinear program to global optimality, here is an example to utilize different sub-solvers with default options set for Alpine:","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"using Alpine\nusing JuMP \nusing Gurobi \nusing Ipopt\n\n# MIP optimizer\nconst gurobi = optimizer_with_attributes(Gurobi.Optimizer, \n                                         MOI.Silent() => true,\n                                         \"Presolve\"   => 1) \n\n# NLP optimizer\nconst ipopt = optimizer_with_attributes(Ipopt.Optimizer, \n                                        MOI.Silent() => true, \n                                        \"sb\" => \"yes\", \n                                        \"max_iter\"   => 9999)\n\n# Global optimizer\nconst alpine = optimizer_with_attributes(Alpine.Optimizer, \n                                         \"nlp_solver\" => ipopt,\n                                         \"mip_solver\" => gurobi)\nm = Model(alpine)","category":"page"},{"location":"choosingsolver/","page":"Choosing Sub-solvers","title":"Choosing Sub-solvers","text":"Similarly, many such optimizer options and nonconvex problems (both NLPs and MINLPs) can be found in the examples folder. ","category":"page"},{"location":"functions/#Functions","page":"Methods","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Methods","title":"Methods","text":"CurrentModule = Alpine","category":"page"},{"location":"functions/#High-level-Algorithmic-Functions","page":"Methods","title":"High-level Algorithmic Functions","text":"","category":"section"},{"location":"functions/","page":"Methods","title":"Methods","text":"These are the high-level algorithmic functions:","category":"page"},{"location":"functions/","page":"Methods","title":"Methods","text":"presolve\nglobal_solve\nlocal_solve\nbounding_solve","category":"page"},{"location":"functions/#Alpine.presolve","page":"Methods","title":"Alpine.presolve","text":"presolve(m::Optimizer)\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.global_solve","page":"Methods","title":"Alpine.global_solve","text":"global_solve(m::Optimizer)\n\nPerform global optimization algorithm that is based on the adaptive piecewise convexification. This iterative algorithm loops over bounding_solve and local_solve until the optimality gap between the lower bound (relaxed problem with min. objective) and the upper bound (feasible problem) is within the user prescribed limits. Each bounding_solve provides a lower bound that serves as the partitioning point for the next iteration (this feature can be modified given a different add_adaptive_partition). Each local_solve provides an incumbent feasible solution. The algorithm terminates when atleast one of these conditions are satisfied: time limit, optimality condition, or iteration limit.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.local_solve","page":"Methods","title":"Alpine.local_solve","text":"Alp.local_solve(m::Optimizer, presolve::Bool=false)\n\nPerform a local NLP or MINLP solve to obtain a feasible solution. The presolve option is set to true when the function is invoked in presolve. Otherwise, the function is invoked from bounding_solve.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.bounding_solve","page":"Methods","title":"Alpine.bounding_solve","text":"Alp.bounding_solve(m::Optimizer; kwargs...)\n\nThis step usually solves a convex MILP/MIQCP/MIQCQP problem for lower bounding the given minimization problem. It solves the problem built upon a piecewise convexification based on the discretization sictionary of some variables. See create_bounding_mip for more details of the problem solved here.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Adapative-Partitioning-Methods","page":"Methods","title":"Adapative Partitioning Methods","text":"","category":"section"},{"location":"functions/","page":"Methods","title":"Methods","text":"create_bounding_mip\npick_disc_vars\nfix_domains\nmin_vertex_cover","category":"page"},{"location":"functions/#Alpine.create_bounding_mip","page":"Methods","title":"Alpine.create_bounding_mip","text":"create_bounding_mip(m::Optimizer; use_disc = nothing)\n\nSet up a MILP bounding model base on variable domain partitioning information stored in use_disc. By default, if use_disc is not provided, it will use m.discretizations store in the Alpine model. The basic idea of this MILP bounding model is to use Tighten McCormick to convexify the original Non-convex region. Among all presented partitions, the bounding model will choose one specific partition as the lower bound solution. The more partitions there are, the better or finer bounding model relax the original MINLP while the more efforts required to solve this MILP is required.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.pick_disc_vars","page":"Methods","title":"Alpine.pick_disc_vars","text":"pickdiscvars(m::Optimizer)\n\nThis function helps pick the variables for discretization. The method chosen depends on user-inputs. In case when indices::Int is provided, the method is chosen as built-in method. Currently, there are two built-in options for users as follows:\n\nmax_cover (Alp.get_option(m, :disc_var_pick)=0, default): pick all variables involved in the non-linear term for discretization\nmin_vertex_cover (Alp.get_option(m, :disc_var_pick)=1): pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified\n\nFor advanced usage, Alp.get_option(m, :disc_var_pick) allows ::Function inputs. User can provide his/her own function to choose the variables for discretization.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.fix_domains","page":"Methods","title":"Alpine.fix_domains","text":"fix_domains(m::Optimizer)\n\nThis function is used to fix variables to certain domains during the local solve process in the global_solve. More specifically, it is used in local_solve to fix binary and integer variables to lower bound solutions and discretizing variables to the active domain according to lower bound solution.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.min_vertex_cover","page":"Methods","title":"Alpine.min_vertex_cover","text":"minvertexcover(m::Optimizer)\n\nmin_vertex_cover chooses the variables based on the minimum vertex cover algorithm for the interaction graph of nonlinear terms which are adaptively partitioned for global optimization. This option can be activated by setting disc_var_pick = 1.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Presolve-Methods","page":"Methods","title":"Presolve Methods","text":"","category":"section"},{"location":"functions/","page":"Methods","title":"Methods","text":"bound_tightening\nminmax_bound_tightening\ncreate_bound_tightening_model\nsolve_bound_tightening_model\nresolve_var_bounds","category":"page"},{"location":"functions/#Alpine.bound_tightening","page":"Methods","title":"Alpine.bound_tightening","text":"bound_tightening(m::Optimizer)\n\nEntry point to the optimization-based bound-tightening (OBBT) algorithm. The aim of the OBBT algorithm is to sequentially tighten the variable bounds until a fixed point is reached.\n\nCurrently, two OBBT methods are implemented in minmax_bound_tightening.\n\n* Bound-tightening with polyhedral relaxations (McCormick, Lambda for convex-hull)\n* Bound-tightening with piecewise polyhedral relaxations: (with three partitions around the local feasible solution)\n\nIf no local feasible solution is obtained, the algorithm defaults to OBBT without partitions\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.minmax_bound_tightening","page":"Methods","title":"Alpine.minmax_bound_tightening","text":"minmax_bound_tightening(m:Optimizer; use_bound::Bool=true, use_tmc::Bool)\n\nThis function implements the OBBT algorithm to tighten the variable bounds. It utilizes either the basic polyhedral relaxations or the piecewise polyhedral relaxations (TMC) to tighten the bounds. The TMC has additional binary variables while performing OBBT.\n\nThe algorithm as two main parameters. The first is the use_tmc, which when set to true invokes the algorithm on the TMC relaxation. The second parameter use_bound takes in the objective value of the local solve solution stored in best_sol for performing OBBT. The use_bound option is set to true when the local solve is successful in obtaining a feasible solution, else this parameter is set to false.\n\nFor details, refer to section 3.1.1 of Nagarajan, Lu, Wang, Bent, Sundar, \"An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs\" link.\n\nSeveral other user-input options can be used to tune the OBBT algorithm. For more details, see Presolve Options.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.create_bound_tightening_model","page":"Methods","title":"Alpine.create_bound_tightening_model","text":"create_bound_tightening_model(m::Optimizer, discretization::Dict, bound::Float64)\n\nThis function takes in the initial discretization information and builds the OBBT model. It is an algorithm specific function called by minmax_bound_tightening.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.solve_bound_tightening_model","page":"Methods","title":"Alpine.solve_bound_tightening_model","text":"solve_bound_tightening_model(m::Optimizer)\n\nA function that solves the min and max OBBT model.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.resolve_var_bounds","page":"Methods","title":"Alpine.resolve_var_bounds","text":"resolve_var_bounds(m::Optimizer)\n\nResolve the bounds of the lifted variable using the information in lvartight and uvartight. This method only takes in known or trivial bounds information to deduce lifted variable bounds and to potentially avoid the cases of infinity bounds.\n\n\n\n\n\nresolve_var_bounds(nonconvex_terms::Dict, discretization::Dict)\n\nFor discretization to be performed, we do not allow a variable being discretized to have infinite bounds.\nThe lifted/auxiliary variables may have infinite bounds and the function infers bounds on these variables. This process\ncan help speed up the subsequent solve times.\n\nOnly used in presolve bound tightening\n\n\n\n\n\n","category":"function"},{"location":"functions/#Utility-Methods","page":"Methods","title":"Utility Methods","text":"","category":"section"},{"location":"functions/","page":"Methods","title":"Methods","text":"update_var_bounds\ndiscretization_to_bounds\ninit_disc\nto_discretization\nflatten_discretization\nadd_adaptive_partition","category":"page"},{"location":"functions/#Alpine.update_var_bounds","page":"Methods","title":"Alpine.update_var_bounds","text":"update_var_bounds(m::Optimizer, discretization::Dict; len::Float64=length(keys(discretization)))\n\nThis function take in a dictionary-based discretization information and convert them into two bounds vectors (lvar, uvar) by picking the smallest and largest numbers. User can specify a certain length that may contains variables that is out of the scope of discretization.\n\nOutput::\n\nl_var::Vector{Float64}, u_var::Vector{Float64}\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.discretization_to_bounds","page":"Methods","title":"Alpine.discretization_to_bounds","text":"discretizationtobounds(d::Dict, l::Int)\n\nSame as update_var_bounds\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.init_disc","page":"Methods","title":"Alpine.init_disc","text":"init_disc(m::Optimizer)\n\nThis function initialize the dynamic discretization used for any bounding models. By default, it takes (.lvarorig, .uvarorig) as the base information. User is allowed to use alternative bounds for initializing the discretization dictionary. The output is a dictionary with MathProgBase variable indices keys attached to the :Optimizer.discretization.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.to_discretization","page":"Methods","title":"Alpine.to_discretization","text":"to_discretization(m::Optimizer, lbs::Vector{Float64}, ubs::Vector{Float64})\n\nUtility functions to convert bounds vectors to Dictionary based structures that is more suitable for partition operations.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.flatten_discretization","page":"Methods","title":"Alpine.flatten_discretization","text":"flatten_discretization(discretization::Dict)\n\nUtility functions to eliminate all partition on discretizing variable and keep the loose bounds.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Alpine.add_adaptive_partition","page":"Methods","title":"Alpine.add_adaptive_partition","text":"add_adaptive_partition(m::Optimizer; use_disc::Dict, use_solution::Vector)\n\nA built-in method used to add a new partition on feasible domains of variables chosen for partitioning.\n\nThis can be illustrated by the following example. Let the previous iteration's partition vector on  variable \"x\" be given by [0, 3, 7, 9]. And say, the lower bounding solution has a value of 4 for variable \"x\". In the case when partition_scaling_factor=4, this function creates the new partition vector as follows: [0, 3, 3.5, 4, 4.5, 7, 9]\n\nThere are two options for this function,\n\n* `use_disc(default=m.discretization)`:: to regulate which is the base to add new partitions on\n* `use_solution(default=m.best_bound_sol)`:: to regulate which solution to use when adding new partitions on\n\nThis function can be accordingly modified by the user to change the behavior of the solver, and thus the convergence.\n\n\n\n\n\n","category":"function"},{"location":"algorithm/#Demo","page":"Demo","title":"Demo","text":"","category":"section"},{"location":"algorithm/","page":"Demo","title":"Demo","text":"CurrentModule = Alpine","category":"page"},{"location":"algorithm/","page":"Demo","title":"Demo","text":"Consider the following problem as an example.","category":"page"},{"location":"algorithm/","page":"Demo","title":"Demo","text":"min_x c^Tx\nst     a_i^Tx text sense_i  b_i forall i\n         l leq x leq u","category":"page"},{"location":"#Alpine's-Documentation","page":"Introduction","title":"Alpine's Documentation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = Alpine","category":"page"},{"location":"#Overview","page":"Introduction","title":"Overview","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Alpine is a Julia/JuMP-based global optimization solver for non-convex programs. Alpine applies a two-stage approach to strengthen piecewise convex relaxations for mixed-integer nonlinear programs (MINLP) with multi-linear terms. In the first stage, Alpine exploits Constraint Programing techniques to contract the variable bounds. Further, it applies optimality-based bound tightening (OBBT) methods iteratively until a fixed point with respect to the bounds are achieved. In the second stage, Alpine partitions the variable domains using an adaptive multivariate partitioning scheme (a parametrized partition at the current solution of sequentially solved MILPs), leading to sparser but tighter relaxations. This approach decouples the number of partitions from the size of the variable domains, leading to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. Alpine also applies polyhedral cutting plane methods to handle convex relaxations of higher-order monomial terms.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Allowable nonlinearities: Alpine can currently handle MINLPs with polynomials in constraints and/or in the objective. Currently, there is no support for exponential cones and Positive Semi-Definite (PSD) cones in MINLPs. Alpine is also a good fit for subsets of the MINLP family, e.g., Mixed-Integer Quadratically Constrainted Quadradic Programs (MIQCQPs), Non-Linear Programs (NLPs), etc.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For more details, check out this video on Alpine at the 2nd Annual JuMP-dev Workshop, held at the Institut de Mathématiques de Bordeaux, June 2018.","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"To use Alpine, first download and install Julia. Note that the current version of Alpine is compatible with Julia 1.0 and later. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The latest stable release of Alpine can be installed using the Julia package manager with","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.add(\"Alpine\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"At least one local nonlinear-programming (NLP) solver and a convex mixed-integer programming (MIP) solver are required for running Alpine. If your problem is an NLP, we recommend Ipopt, and else if your problem is a mixed-integer NLP (MINLP), we recommend Juniper. For the underlying MIP solver, Gurobi is highly recommended, as it is fast, scaleable and can be used to solve on fairly large-scale non-convex programs. However, the open-source HiGHS or GLPK solver is also compatible with Alpine. For example, Ipopt and Gurobi can be installed via the package manager with","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.add(\"Ipopt\")\nPkg.add(\"Gurobi\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For more details on choosing sub-solvers for Alpine, check out here. ","category":"page"},{"location":"#Unit-Tests","page":"Introduction","title":"Unit Tests","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"To run the tests in the package, run the following command after installing the Alpine package.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.test(\"Alpine\")","category":"page"},{"location":"#Citing-Alpine","page":"Introduction","title":"Citing Alpine","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"If you find Alpine useful in your work, we request you to cite the following paper [pdf] [pdf]: ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"@article{alpine_JOGO2019,\n  title = {An adaptive, multivariate partitioning algorithm for global optimization of nonconvex programs},\n  author = {Nagarajan, Harsha and Lu, Mowen and Wang, Site and Bent, Russell and Sundar, Kaarthik},\n  journal = {Journal of Global Optimization},\n  year = {2019},\n  issn = {1573-2916},\n  doi = {10.1007/s10898-018-00734-1},\n}\n\n@inproceedings{alpine_CP2016,\n  title = {Tightening {McCormick} relaxations for nonlinear programs via dynamic multivariate partitioning},\n  author = {Nagarajan, Harsha and Lu, Mowen and Yamangil, Emre and Bent, Russell},\n  booktitle = {International Conference on Principles and Practice of Constraint Programming},\n  pages = {369--387},\n  year = {2016},\n  organization = {Springer},\n  doi = {10.1007/978-3-319-44953-1_24},\n}","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"If you find the underlying piecewise polyhedral formulations implemented in Alpine useful in your work, we kindly request that you cite the following papers (link-1, link-2): ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"@article{alpine_ORL2021,\n  title = {Piecewise polyhedral formulations for a multilinear term},\n  author = {Sundar, Kaarthik and Nagarajan, Harsha and Linderoth, Jeff and Wang, Site and Bent, Russell},\n  journal = {Operations Research Letters},\n  volume = {49},\n  number = {1},\n  pages = {144--149},\n  year = {2021},\n  publisher = {Elsevier}\n}\n\n@article{alpine_OptOnline2022,\n    title={Piecewise Polyhedral Relaxations of Multilinear Optimization},\n    author={Kim, Jongeun and Richard, Jean-Philippe P. and Tawarmalani, Mohit},\n    eprinttype={Optimization Online},\n    date={2022}\n}","category":"page"},{"location":"parameters/#Solver-options-for-Alpine","page":"Solver Options","title":"Solver options for Alpine","text":"","category":"section"},{"location":"parameters/","page":"Solver Options","title":"Solver Options","text":"CurrentModule = Alpine","category":"page"},{"location":"parameters/#General-Options","page":"Solver Options","title":"General Options","text":"","category":"section"},{"location":"parameters/","page":"Solver Options","title":"Solver Options","text":"Here are a few general solver options which control the performance of Alpine:","category":"page"},{"location":"parameters/","page":"Solver Options","title":"Solver Options","text":"log_level (default = 0): verbosity level of Alpine; choose 1 for turning on logging, else 100 for detailed debugging mode.\ntime_limit (default = Inf): total limit on run time for Alpine in seconds.\nmax_iter (default = 999): total number of iterations allowed in the global_solve. \nrel_gap (default = 1e-4): relative gap considered for global convergence during global_solve. Bounds are evaluated using fracUB-LBUB cdot 100 .\ntol (default = 1e-6): numerical tolerance used during the process of global_solve. ","category":"page"},{"location":"parameters/#Adaptive-Partitioning-Options","page":"Solver Options","title":"Adaptive Partitioning Options","text":"","category":"section"},{"location":"parameters/","page":"Solver Options","title":"Solver Options","text":"apply_partitioning (default = true): applies Alpine's built-in MIP-based partitioning algorithm only when activated; else terminates with the presolve solution.\npartition_scaling_factor (default = 10): used during add_adaptive_partition for scaling the width of new partitions relative to the active partition chosen in the sequentially-solved lower-bounding MIP models. This value can substantially affect the run time for global convergence; this value can be set to different integer values (>= 4) for various classes of problems. \ndisc_var_pick (default = 0): controls Alpine's algorithm used for selecting variables for partitioning; 0 is for max-cover, 1 is for minimum-vertex-cover. This parameter allows functional inputs.\ndisc_add_partition_method: allows functional input on how new partitions could be constructed.","category":"page"},{"location":"parameters/#Presolve-Options","page":"Solver Options","title":"Presolve Options","text":"","category":"section"},{"location":"parameters/","page":"Solver Options","title":"Solver Options","text":"presolve_track_time (default = true): includes/excludes presolve run time in the total Alpine's run time. \npresolve_bt (default = false): performs sequential, optimization-based bound tightening (OBBT) at the presolve step. \npresolve_bt_max_iter (default = 9999): maximum number of iterations allowed using the sequential OBBT step. \npresolve_bt_width_tol (default = 1e-3): numerical tolerance value used in the OBBT step. Note that smaller values of this tolerance can lead to inaccurate solutions.  \npresolve_bt_algo (default = 1): method chosen to perform the presolve step; choose 1 for the built-in OBBT, else 2 for user-input functions for presolve.\npresolve_bt_relax_integrality (default = false): relaxes the integrality of the existing integer variables in the OBBT step, if the input problem is an MINLP. \npresolve_bt_mip_time_limit (default = Inf): time limit for individual MILPs solved during the sequential OBBT procedure. ","category":"page"},{"location":"parameters/","page":"Solver Options","title":"Solver Options","text":"Note that the above-mentioned list of solver options is not comprehensive, but can be found in solver.jl. ","category":"page"}]
}
