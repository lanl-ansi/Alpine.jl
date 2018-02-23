"""
    Collects distance of each variable's bounding solution to best feasible solution and select the ones that is the furthest
    Currently don't support recursively convexification
"""
function update_disc_cont_var(m::PODNonlinearModel)

    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    # If no feasible solution found, do NOT update
    if m.status[:feasible_solution] != :Detected
        println("no feasible solution detected. No update disc var selection.")
        return
    end

	var_idxs = copy(m.candidate_disc_vars)
    var_diffs = Vector{Float64}(m.num_var_orig+length(keys(m.linear_terms))+length(keys(m.nonconvex_terms)))

    for i in 1:m.num_var_orig       # Original Variables
        var_diffs[i] = abs(m.best_sol[i]-m.best_bound_sol[i])
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
    weighted_min_vertex_cover(m, distance)

    (m.loglevel > 100) && println("updated partition var selection => $(m.disc_vars)")
    return
end

function update_disc_int_var(m::PODNonlinearModel)

    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    # If all discretized integer solution have been fully discovered, then include all other integer variables
    checker = [is_fully_discovered_integer(i, m.discretization[i]) for i in m.disc_vars if m.var_type[i] == :Int]
    if prod(checker)
        for i in m.candidate_disc_vars
            if m.var_type[i] == :Int && !(i in m.disc_vars)
                m.loglevel > 99 && println("PUMPing Integer VAR$(i) into discretization var set.")
                push!(m.disc_vars, i)
            end
        end
    end

    return
end

"""
    One-time rounding heuristic to obtain a feasible solution
"""
function heu_basic_rounding(m::PODNonlinearModel, local_model)

    println("Basic Rounding Heuristic Avtivated...")

    convertor = Dict(:Max=>:>, :Min=>:<)

    rounded_sol = round_sol(m, local_model)
    l_var, u_var = fix_domains(m, discrete_sol = rounded_sol)

    heuristic_model = interface_init_nonlinear_model(m.nlp_solver)
    interface_load_nonlinear_model(m, heuristic_model, l_var, u_var)
    interface_optimize(heuristic_model)
    heuristic_model_status = interface_get_status(heuristic_model)

    if heuristic_model_status in [:Infeasible, :Error]
        println("Rounding obtained an Infeasible point.")
        push!(m.logs[:obj], "-")
        return :Infeasibles
    elseif heuristic_model_status in [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
        candidate_obj = interface_get_objval(heuristic_model)
        candidate_sol = round.(interface_get_solution(heuristic_model), 5)
        update_incumb_objective(m, candidate_obj, candidate_sol)
        println("Rounding obtained a feasible solution OBJ = $(m.best_obj)")
        return :LocalOptimal
    else
        error("[EXCEPTION] Unknown NLP solver status.")
    end

    return
end

function heu_pool_boosting(m::PODNonlinearModel, local_model)

    println("Heuristic Local Search : LB pool-based multi-search")

    convertor = Dict(:Max=>:>, :Min=>:<)


    return
end
