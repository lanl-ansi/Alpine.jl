var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#POD.jl-Documentation-1",
    "page": "Home",
    "title": "POD.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "POD.jl relies on piecewise convex relaxations to solve mixed-integer nonlinear programs (MINLP). The main algorithmic idea is to strenghten the convex relaxation by systematically partitioning variable domain and use the relaxation to guide the local search. Instead of equally partitioning the domains of variables appearing in multi-linear terms (predominantly common in the literature), we construct sparser partitions yet tighter relaxations by adaptively partitioning the variable domains in regions of interest. This approach decouples the number of partitions from the size of the variable domains, leads to a significant reduction in computation time, and limits the number of binary variables that are introduced by the partitioning. To describe the piecewise function, POD.jl deploy a locally idea convexification formulation with performance gurantee in practical computation. Furthermore, POD.jl tries constraint programing techniques to contract the variable bounds for further performance improvements. POD.jl is specifically designed to fit Julia JuMP model. It performs expression-based analyze to structurally understand the problem and apply the algorithm smartly. The code is design to embrace the open-source community with API-based inputs, which is used to handle problems of users' interests with user's idea bridged with POD's algorithm."
},

{
    "location": "index.html#Documentation-Structure-1",
    "page": "Home",
    "title": "Documentation Structure",
    "category": "section",
    "text": "We kindly ask the user to walk through the \"Getting Started\" section for a breif, basic usage guidance on Installation, Usage, and Choosing Sub-Solvers. The list of Parameters is documented and it can be used to control the solver. For more advanced control of the solver, user should consult with Expression Guideline and Hacking POD for more details about how different expressions can affect the solver behavior and how to apply user APIs."
},

{
    "location": "index.html#Important-Note-1",
    "page": "Home",
    "title": "Important Note",
    "category": "section",
    "text": "POD is a developing software as well as a developing research project. POD is not yet numerically robust and may require dedicated attention when encountering new instances. Hence, POD is designed to accept user APIs for algorithm tuning. If you need such attention, please let us know."
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
    "text": "Currently, POD.jl is a private repository that is owned by LANL-ANSI group. To install, doPkg.clone(\"https://github.com/lanl-ansi/POD.jl.git\")To developers: it is highly recommend that any further development on POD is conducted on a forked repo."
},

{
    "location": "usage.html#",
    "page": "Basic Usage",
    "title": "Basic Usage",
    "category": "page",
    "text": ""
},

{
    "location": "usage.html#Usage-1",
    "page": "Basic Usage",
    "title": "Usage",
    "category": "section",
    "text": ""
},

{
    "location": "usage.html#Building-a-problem-1",
    "page": "Basic Usage",
    "title": "Building a problem",
    "category": "section",
    "text": "For simplicity, we present a simple nonlinear program without any discrete variable for demonstration. The model is implemented in JuMP format. Two opon-source sub-solvers, Ipopt and Cbc, are used to handle local NLP (nlp_local_solver) and relaxed MIP (mip_solver). We control POD.jl to provide basic solving log through log_level. Before you solve, make sure that the two sub-solvers are properly installed and built.using JuMP, POD\nusing Ipopt, Cbc\n\nm = Model(solver=PODSolver(nlp_local_solver = IpoptSolver(print_level=0),\n                           mip_solver = CbcSolver(logLevel=0),\n                           rel_gap = 0.001,\n                           log_level = 1))\n\n@variable(m, x[1:8])\n\nsetlowerbound(x[1], 100)\nsetlowerbound(x[2], 1000)\nsetlowerbound(x[3], 1000)\nsetlowerbound(x[4], 10)\nsetlowerbound(x[5], 10)\nsetlowerbound(x[6], 10)\nsetlowerbound(x[7], 10)\nsetlowerbound(x[8], 10)\n\nsetupperbound(x[1], 10000)\nsetupperbound(x[2], 10000)\nsetupperbound(x[3], 10000)\nsetupperbound(x[4], 1000)\nsetupperbound(x[5], 1000)\nsetupperbound(x[6], 1000)\nsetupperbound(x[7], 1000)\nsetupperbound(x[8], 1000)\n\n@constraint(m, 0.0025*(x[4]+x[6]) <= 1)\n@constraint(m, 0.0025*(x[5]-x[4]+x[7]) <= 1)\n@constraint(m, 0.01(x[8]-x[5]) <= 1)\n@NLconstraint(m, 100*x[1] - x[1]*x[6] + 833.33252*x[4] <= 83333.333)\n@NLconstraint(m, x[2]*x[4] - x[2]*x[7] - 1250*x[4] + 1250*x[5] <= 0)\n@NLconstraint(m, x[3]*x[5] - x[3]*x[8] - 2500*x[5] + 1250000 <= 0)\n\n@objective(m, Min, x[1]+x[2]+x[3])\n\nsolve(m)"
},

{
    "location": "usage.html#Interpreting-POD-log-1",
    "page": "Basic Usage",
    "title": "Interpreting POD log",
    "category": "section",
    "text": "POD.jl will output the log below. The output contains two main sections: 1) basic problem dimension and solver configuration summaries, and 2)presolve and main algorithms process. In the latter, The first two columns NLP and MIP are local solves and realxation solve objective value. Column Objective and Bound tracks the incumbent upper and lower bound of the problem. Without the loss of generality, upper bound is always defined through a feasible solution. Column GAP% is relative optimality gap evaluated through,    textbfGap = fracUB-LB+UB times 100Total wall time performance is tracked through column CLOCK. Last but not least, the iteration is tracked incrementally in column Iter.full problem loaded into POD\nproblen sense Min\nnumber of constraints = 6\nnumber of non-linear constraints = 3\nnumber of linear constraints = 3\nnumber of variables = 8\nMINLP local solver = POD\nNLP local solver = Ipopt\nMIP solver = Cbc\nmaximum solution time = Inf\nmaximum iterations =  99\nrelative optimality gap criteria = 0.00100 (0.1000 %)\ndetected nonlinear terms = 5\nnumber of variables involved in nonlinear terms = 8\nnumber of selected variables to discretize = 8\nbilinear treatment = convex hull formulation\nmonomial treatment = convex hull formulation\nusing convex hull : sos2 formulation\nusing method 0 for picking discretization variable...\nadaptively adding discretization ratio = 4\n\nPOD algorithm presolver started.\nPerforming local solve to obtain a feasible solution.\nPresolve ended.\nPresolve time = 0.02s\nm.logs[:time_left] = Inf\n | NLP           | MIP           || Objective     | Bound         | GAP%          | CLOCK         | | Iter\n | 7049.2479     | 3004.247      || 7049.2479     | 3004.247      | 57.38202      | 0.07s         | | 1\n | -             | 4896.6075     || 7049.2479     | 4896.6075     | 30.53716      | 0.33s         | | 2\n | 7049.2479     | 5871.5307     || 7049.2479     | 5871.5307     | 16.70699      | 0.75s         | | 3\n | -             | 6717.2923     || 7049.2479     | 6717.2923     | 4.70909       | 2.8s          | | 4\n | 7049.267      | 6901.8131     || 7049.2479     | 6901.8131     | 2.0915        | 8.55s         | | 5\n | 7049.9758     | 7020.723      || 7049.2479     | 7020.723      | 0.40465       | 15.81s        | | 6\n | 7050.6454     | 7042.4094     || 7049.2479     | 7042.4094     | 0.09701       | 31.68s        | | 7\n | 7050.7924     | 7042.4094     || 7049.2479     | 7042.4094     | 0.09701       | 31.69s        | | finish\n\n POD ended with status Optimal\n:Optimal"
},

{
    "location": "choosingsolver.html#",
    "page": "Choosing Sub-Solvers",
    "title": "Choosing Sub-Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "choosingsolver.html#Choosing-Sub-Solvers-1",
    "page": "Choosing Sub-Solvers",
    "title": "Choosing Sub-Solvers",
    "category": "section",
    "text": "The design of the POD.jl requires a variety of programming problems to be solved underneath the surface. More importantly, the potential of POD.jl is best exploited when utilizing strong sub-solvers. The performance of the example shown in Usage is approximatly 1/10 if MIP solver is Gurobi rather than Cbc. POD.jl will customize the algorithm when utilizing different sub-solver. For example, POD.jl will utilize tune starting values for MIP sub-solver and resolve infeasibility with local NLP/MINLP solvers. For best algorithmic performance, it's recommend that dedicated solvers to be used for these operations. Currently, the following sub-solvers with Julia Interface is partially supported by POD.jl:Solver Julia Package Target Program\nCPLEX CPLEX.jl MILP, MIQCP\nCbc Cbc.jl MILP\nGurobi Gurobi.jl MILP, MIQCP\nIpopt Ipopt.jl Local NLP\nBonmin Bonmin.jl Local MINLP\nArtelys KNITRO KNITRO.jl Local MINLPAs the development of AMP contineous, supports fo Mosek, GLPK, NLopt, Xpress is already scheduled for the roadmap of POD.jl."
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
    "text": ""
},

{
    "location": "parameters.html#General-Solver-Control-1",
    "page": "Parameters",
    "title": "General Solver Control",
    "category": "section",
    "text": "Keyword Type Description\nlog_level Int ...\ntimeout Float ...\nmax_iter Int ...\nrel_gap Float ...\ntol Float ..."
},

{
    "location": "parameters.html#Sub-solvers-Control-1",
    "page": "Parameters",
    "title": "Sub-solvers Control",
    "category": "section",
    "text": "Keyword Type Description\nminlp_local_solver AbstractMathProgSolver ...\nnlp_local_solver AbstractMathProgSolver ...\nmip_solver AbstractMathProgSolver ..."
},

{
    "location": "parameters.html#Main-Algorithm-Control-1",
    "page": "Parameters",
    "title": "Main Algorithm Control",
    "category": "section",
    "text": "Keyword Type Description\npresolve Int ...\npresolve_algo Int ...\npresolve_max_iter Int ...\npresolve_mip_relax Int ...\npresolve_timeout Int ...\nrecognize_convex Bool ...\ndisc_ratio Int ...\ndisc_var_pick_algo Int ...\ndisc_abs_width_tol Int ...\ndisc_consec_forbid Int ...\nbounds_propagation Bool ..."
},

{
    "location": "parameters.html#User-API-Control-1",
    "page": "Parameters",
    "title": "User API Control",
    "category": "section",
    "text": "Keyword Type Description\nnlterm_patterns Any ...\nconvex_patterns Any ...\nconvexify_methods Any ...\ndisc_partition_method Any ..."
},

{
    "location": "expression.html#",
    "page": "Expression Guideline",
    "title": "Expression Guideline",
    "category": "page",
    "text": ""
},

{
    "location": "expression.html#Expression-Guideline-1",
    "page": "Expression Guideline",
    "title": "Expression Guideline",
    "category": "section",
    "text": ""
},

{
    "location": "hacking.html#",
    "page": "Hacking POD",
    "title": "Hacking POD",
    "category": "page",
    "text": ""
},

{
    "location": "hacking.html#Hacking-POD-1",
    "page": "Hacking POD",
    "title": "Hacking POD",
    "category": "section",
    "text": "This page give more detailed tutorial on how to hack this solver for different behaviors..."
},

]}
