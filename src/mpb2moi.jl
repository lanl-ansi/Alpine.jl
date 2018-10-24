# This file carries the goal of getting ready to transit from MPB interface to MOI interface
# We wrap the functions
function interface_init_nonlinear_model(solver::Any)
    return MathProgBase.NonlinearModel(solver)
end

function interface_init_nonlinear_data(nonlinear_model::MathProgBase.AbstractNLPEvaluator)
    MathProgBase.initialize(nonlinear_model, [:Grad, :Jac, :Hess, :HessVec, :ExprGraph]) # Safety scheme for sub-solvers re-initializing the NLPEvaluator
    return
end

function interface_load_nonlinear_model(m::PODNonlinearModel, nonlinear_model::Any, l_var::Vector, u_var::Vector)
    MathProgBase.loadproblem!(nonlinear_model, m.num_var_orig,
                                               m.num_constr_orig,
                                               l_var,
                                               u_var,
                                               m.l_constr_orig,
                                               m.u_constr_orig,
                                               m.sense_orig,
                                               m.d_orig)

    return
end

function interface_set_warmstart(nonlinear_model::Any, startval::Vector)
    MathProgBase.setwarmstart!(nonlinear_model, startval)
    return
end

function interface_set_vartype(nonlinear_model::Any, vartype::Vector)
    MathProgBase.setvartype!(nonlinear_model, vartype)
    return
end

function interface_optimize(nonlinear_model::Any)
    MathProgBase.optimize!(nonlinear_model)
    return
end

function interface_get_status(nonlinear_model::Any)
    return MathProgBase.status(nonlinear_model)
end

function interface_get_solution(nonlinear_model::Any)
    return MathProgBase.getsolution(nonlinear_model)
end

function interface_get_objval(nonlinear_model::Any)
    return MathProgBase.getobjval(nonlinear_model)
end

function interface_get_obj_expr(nonlinear_model::MathProgBase.AbstractNLPEvaluator)
    return MathProgBase.obj_expr(nonlinear_model)
end

function interface_get_constr_expr(nonlinear_model::MathProgBase.AbstractNLPEvaluator, idx::Int)
    return MathProgBase.constr_expr(nonlinear_model, idx)
end

function interface_is_constr_linear(nonlinear_model::MathProgBase.AbstractNLPEvaluator, idx::Int)
    return MathProgBase.isconstrlinear(nonlinear_model, idx)
end

function interface_is_obj_linear(nonlinear_model::MathProgBase.AbstractNLPEvaluator)
    return MathProgBase.isobjlinear(nonlinear_model)
end

function interface_eval_g(nonlinear_model::MathProgBase.AbstractNLPEvaluator, rhs, sol)
    MathProgBase.eval_g(nonlinear_model, rhs, sol)
    return
end

function interface_get_rawsolver(m::JuMP.Model)
    return MathProgBase.getrawsolver(internalmodel(m))
end

MathProgBase.setwarmstart!(m::PODNonlinearModel, x) = (m.var_start_orig = x)
MathProgBase.setvartype!(m::PODNonlinearModel, v::Vector{Symbol}) = (m.var_type_orig = v)
MathProgBase.status(m::PODNonlinearModel) = m.pod_status
MathProgBase.getobjval(m::PODNonlinearModel) = m.best_obj
MathProgBase.getobjbound(m::PODNonlinearModel) = m.best_bound
MathProgBase.getsolution(m::PODNonlinearModel) = m.best_sol
MathProgBase.getsolvetime(m::PODNonlinearModel) = m.logs[:total_time]
