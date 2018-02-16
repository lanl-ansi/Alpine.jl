"""
    Evaluate a solution feasibility: Solution bust be in the feasible category and evaluated rhs must be feasible
"""
function eval_feasibility(m::PODNonlinearModel, sol::Vector)

    length(sol) == m.num_var_orig || error("Candidate solution length mismatch.")


    for i in 1:m.num_var_orig
        # Check solution category and bounds
        if m.var_type[i] == :Bin
            isapprox(sol[i], 1.0;atol=m.tol) || isapprox(sol[i], 0.0;atol=m.tol) || return false
        elseif m.var_type[i] == :Int
            isapprox(mod(sol[i], 1.0), 0.0;atol=m.tol) || return false
        end
        # Check solution bounds (with tight bounds)
        sol[i] <= m.l_var_tight[i] - m.tol || return false
        sol[i] >= m.u_var_tight[i] + m.tol || return false
    end

    # Check constraint violation
    eval_rhs = zeros(m.num_constr_orig)
    MathProgBase.eval_g(m.d_orig, eval_rhs, rounded_sol)
    feasible = true
    for i in 1:m.num_constr_orig
        if m.constr_type_orig[i] == :(==)
            if !isapprox(eval_rhs[i], m.l_constr_orig[i]; atol=m.tol)
                feasible = false
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) != RHS $(m.l_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        elseif m.constr_type_orig[i] == :(>=)
            if !(eval_rhs[i] >= m.l_constr_orig[i] - m.tol)
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !>= RHS $(m.l_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        elseif m.constr_type_orig[i] == :(<=)
            if !(eval_rhs[i] <= m.u_constr_orig[i] + m.tol)
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !<= RHS $(m.u_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        end
    end

    return feasible
end

function round_sol(m::PODNonlinearModel, nlp_model)
    relaxed_sol = MathProgBase.getsolution(nlp_model)
    return [m.var_type_orig[i] in [:Bin, :Int] ? round(relaxed_sol[i]) : relaxed_sol[i] for i in 1:m.num_var_orig]
end
