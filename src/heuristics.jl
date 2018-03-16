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
    For integer solutions
"""
function heu_basic_rounding(m::PODNonlinearModel, local_model)

    println("Basic Rounding Heuristic Avtivated...")

    convertor = Dict(:Max=>:>, :Min=>:<)

    rounded_sol = round_sol(m, nlp_model=local_model)
    l_var, u_var = fix_domains(m, discrete_sol = rounded_sol)

    heuristic_model = interface_init_nonlinear_model(m.nlp_solver)
    interface_load_nonlinear_model(m, heuristic_model, l_var, u_var)
    interface_optimize(heuristic_model)
    heuristic_model_status = interface_get_status(heuristic_model)

    if heuristic_model_status in [:Infeasible, :Error]
        m.loglevel > 0 && println("Rounding obtained an Infeasible point.")
        push!(m.logs[:obj], "INF")
        return :Infeasibles
    elseif heuristic_model_status in [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
        candidate_obj = interface_get_objval(heuristic_model)
        candidate_sol = round.(interface_get_solution(heuristic_model), 5)
        update_incumb_objective(m, candidate_obj, candidate_sol)
        m.loglevel > 0 && println("Rounding obtained a feasible solution OBJ = $(m.best_obj)")
        return :LocalOptimal
    else
        error("[EXCEPTION] Unknown NLP solver status.")
    end

    return
end

"""
    Use all lower bound solution pools as starting points
"""
function heu_pool_multistart(m::PODNonlinearModel)

    convertor = Dict(:Max=>:>, :Min=>:<)
    m.sense_orig == :Min ? incumb_obj = Inf : incumb_obj = -Inf
    incumb_sol = []
    found_feasible = false

    for i in 1:m.bound_sol_pool[:cnt]
        if !m.bound_sol_pool[:ubstart][i]
            rounded_sol = round_sol(m, nlp_sol=m.bound_sol_pool[:sol][i])
            l_var, u_var = fix_domains(m, discrete_sol=rounded_sol, use_orig=true)
            heuristic_model = interface_init_nonlinear_model(m.nlp_solver)
            interface_load_nonlinear_model(m, heuristic_model, l_var, u_var)
            interface_optimize(heuristic_model)
            heuristic_model_status = interface_get_status(heuristic_model)
            if heuristic_model_status in [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
                candidate_obj = interface_get_objval(heuristic_model)
                if eval(convertor[m.sense_orig])(candidate_obj, incumb_obj)
                    incumb_obj = candidate_obj
                    incumb_sol = round.(interface_get_solution(heuristic_model), 5)
                    m.loglevel > 0 && println("Feasible solution obtained using lower bound solution pool [SOL:$(i)] [OBJ=$(incumb_obj)]")
                end
                found_feasible = true
            else
                m.loglevel > 99 && println("Multi-start heuristic returns $(heuristic_model_status) [SOL:$(i)]")
            end
            m.bound_sol_pool[:ubstart][i] = true
        end
    end

    if found_feasible
        update_incumb_objective(m, incumb_obj, incumb_sol)
        return :LocalOptimal
    end

    return :Infeasible
end
