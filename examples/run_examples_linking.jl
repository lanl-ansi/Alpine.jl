# - Necessary - #
using Alpine
using JuMP
using Gurobi
using Ipopt
using Juniper
using Logging

# Install library package. To install, at your Julia command prompt,
#    Pkg.add(path="https://github.com/lanl-ansi/MINLPLib.jl.git")
# Reference: https://github.com/lanl-ansi/MINLPLib.jl
import MINLPLib

# Choose underlying solvers for Alpine
ipopt = get_ipopt()
gurobi = get_gurobi()
juniper = get_juniper(gurobi, ipopt)

#= Global solver
 Hints: 
 => Try different integer values (>=4) for `disc_ratio` to potentially observe 
    better Alpine run times (can be instance specific).
 => Choose `presolve_bt` to `false` if you prefer the bound tightening (OBBT) presolve to be turned off. 
 => If you prefer to use Alpine for only OBBT presolve, without any paritioning applied to the 
    nonlinear terms, include option "apply_partitioning" below and set it to false. 
=#

function run(use_linking_constraints, time_limit, instance_name)
    alpine = JuMP.optimizer_with_attributes(
        Alpine.Optimizer,
        "minlp_solver" => juniper,
        "nlp_solver" => ipopt,
        "mip_solver" => gurobi,
        "presolve_bt" => true,
        "disc_ratio" => 10,
        "time_limit" => time_limit,
        "linking_constraints" => use_linking_constraints,
    )
    m = MINLPLib.fetch_model(instance_name)
    JuMP.set_optimizer(m, alpine)

    t_begin = time()
    JuMP.optimize!(m)
    objval = JuMP.objective_value(m)
    soltime = time() - t_begin

    return objval, soltime
end

# Linking constraints are only added when there are multilinear terms whose degree is at least 3.
# Choose instances in `mult3`` or `mult4`` libraries to see the improvement. 
const instance_name = "mult3/m_10_3_10_100_1"
const time_limit = 60

objval_nolink, soltime_nolink = run(false, time_limit, instance_name)
objval_link, soltime_link = run(true, time_limit, instance_name)

@info "Solution Time:\t link = $soltime_link, nolink = $soltime_nolink"
@info "Objective Value:\t link = $objval_link, nolink = $objval_nolink"
