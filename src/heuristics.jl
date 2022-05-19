"""
    Ranking of variables involved in nonlinear terms for piecewise adaptive partitioning:
    Ranked based on the absolute gap between each variable's solution from the lower-bounding MIP and the best feasible solution to the MINLP.
    Currently doesn't support recursive convexification
"""
function update_disc_cont_var(m::Optimizer)

    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    # If no feasible solution is found, do NOT update
    if !m.detected_feasible_solution
        println("no feasible solution detected. No update disc var selection.")
        return
    end

	var_idxs = copy(m.candidate_disc_vars)
    var_diffs = Vector{Float64}(undef, m.num_var_orig+length(keys(m.linear_terms))+length(keys(m.nonconvex_terms)))

    for i in 1:m.num_var_orig       # Original Variables
        var_diffs[i] = abs(m.best_sol[i] - m.best_bound_sol[i])
    end

    for i in 1:length(keys(m.linear_terms))
        for j in keys(m.linear_terms)   # sequential evaluation to avoid dependency issue
            if m.linear_terms[j][:id] == i
                var_diffs[m.linear_terms[j][:lifted_var_ref].args[2]] = m.linear_terms[j][:evaluator](m.linear_terms[j], var_diffs)
            end
        end
    end

    for i in 1:length(keys(m.nonconvex_terms))
        for j in keys(m.nonconvex_terms)    # sequential evaluation to avoid dependency issue
            if m.nonconvex_terms[j][:id] == i
                var_diffs[m.nonconvex_terms[j][:lifted_var_ref].args[2]] = m.nonconvex_terms[j][:evaluator](m.nonconvex_terms[j], var_diffs)
            end
        end
    end

    distance = Dict(zip(var_idxs,var_diffs))
    Alp.weighted_min_vertex_cover(m, distance)

    (Alp.get_option(m, :log_level) > 100) && println("updated partition var selection => $(m.disc_vars)")
    return
end

function update_disc_int_var(m::Optimizer)

    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    # Additional error checking scheme
    :Int in m.var_type && error("Alpine does not support MINLPs with geenric integer (non-binary) variables yet! Try Juniper.jl for finding a local feasible solution")

    return
end

"""
    Use solutions from the MIP solution pool as starting points
"""
function heu_pool_multistart(m::Optimizer)

    convertor = Dict(MOI.MAX_SENSE => :>, MOI.MIN_SENSE => :<)
    Alp.is_min_sense(m) ? incumb_obj = Inf : incumb_obj = -Inf
    incumb_sol = []
    found_feasible = false

    for i in 1:m.bound_sol_pool[:cnt]
        if !m.bound_sol_pool[:ubstart][i]

            rounded_sol = Alp.round_sol(m, m.bound_sol_pool[:sol][i])
            l_var, u_var = Alp.fix_domains(m, discrete_sol=rounded_sol, use_orig=true)
            heuristic_model = MOI.instantiate(Alp.get_option(m, :nlp_solver), with_bridge_type=Float64)
            x = Alp.load_nonlinear_model(m, heuristic_model, l_var, u_var)
            
            MOI.optimize!(heuristic_model)

            heuristic_model_status = MOI.get(heuristic_model, MOI.TerminationStatus())
            if heuristic_model_status in STATUS_OPT || heuristic_model_status in STATUS_LIMIT
                candidate_obj = MOI.get(heuristic_model, MOI.ObjectiveValue())
                if eval(convertor[m.sense_orig])(candidate_obj, incumb_obj)
                    incumb_obj = candidate_obj
                    incumb_sol = round.(MOI.get(heuristic_model, MOI.VariablePrimal(), x), 5)
                    Alp.get_option(m, :log_level) > 0 && println("Feasible solution obtained using lower bound solution pool [SOL:$(i)] [OBJ=$(incumb_obj)]")
                end
                found_feasible = true
            else
                Alp.get_option(m, :log_level) > 99 && println("Multi-start heuristic returns $(heuristic_model_status) [SOL:$(i)]")
            end
            m.bound_sol_pool[:ubstart][i] = true
        end
    end

    if found_feasible
        Alp.update_incumb_objective(m, incumb_obj, incumb_sol)
        return MOI.LOCALLY_SOLVED
    end

    return MOI.LOCALLY_INFEASIBLE
end


#-----------------------------------------------------------------#
#                    UNSUPPORTED FUNCTIONS                        #
#-----------------------------------------------------------------#

# """
#     One-time rounding heuristic to obtain a feasible solution
#     For integer solutions
# """
# function heu_basic_rounding(m::Optimizer, relaxed_sol)

#     println("Basic Rounding Heuristic Activated...")

#     convertor = Dict(MOI.MAX_SENSE => :>, MOI.MIN_SENSE => :<)

#     rounded_sol = round_sol(m, relaxed_sol)
#     l_var, u_var = fix_domains(m, discrete_sol = rounded_sol)

#     heuristic_model = MOI.instantiate(Alp.get_option(m, :nlp_solver), with_bridge_type=Float64)
#     x = load_nonlinear_model(m, heuristic_model, l_var, u_var)
#     MOI.optimize!(heuristic_model)
#     heuristic_model_status = MOI.get(heuristic_model, MOI.TerminationStatus())

#     if heuristic_model_status == MOI.OTHER_ERROR || heuristic_model_status in STATUS_INF
#         Alp.get_option(m, :log_level) > 0 && println("Rounding obtained an Infeasible point.")
#         push!(m.logs[:obj], "INF")
#         return MOI.LOCALLY_INFEASIBLE
#     elseif heuristic_model_status in STATUS_OPT || heuristic_model_status in STATUS_LIMIT
#         candidate_obj = MOI.get(heuristic_model, MOI.ObjectiveValue())
#         candidate_sol = round.(MOI.get(heuristic_model, MOI.VariablePrimal(), x), 5)
#         update_incumb_objective(m, candidate_obj, candidate_sol)
#         Alp.get_option(m, :log_level) > 0 && println("Rounding obtained a feasible solution OBJ = $(m.best_obj)")
#         return MOI.LOCALLY_SOLVED
#     else
#         error("[EXCEPTION] Unknown NLP solver status.")
#     end

#     return
# end