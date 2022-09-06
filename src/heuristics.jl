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
    var_diffs = Vector{Float64}(
        undef,
        m.num_var_orig + length(keys(m.linear_terms)) + length(keys(m.nonconvex_terms)),
    )

    for i in 1:m.num_var_orig       # Original Variables
        var_diffs[i] = abs(m.best_sol[i] - m.best_bound_sol[i])
    end

    for i in 1:length(keys(m.linear_terms))
        for j in keys(m.linear_terms)   # sequential evaluation to avoid dependency issue
            if m.linear_terms[j][:id] == i
                var_diffs[m.linear_terms[j][:lifted_var_ref].args[2]] =
                    m.linear_terms[j][:evaluator](m.linear_terms[j], var_diffs)
            end
        end
    end

    for i in 1:length(keys(m.nonconvex_terms))
        for j in keys(m.nonconvex_terms)    # sequential evaluation to avoid dependency issue
            if m.nonconvex_terms[j][:id] == i
                var_diffs[m.nonconvex_terms[j][:lifted_var_ref].args[2]] =
                    m.nonconvex_terms[j][:evaluator](m.nonconvex_terms[j], var_diffs)
            end
        end
    end

    distance = Dict(zip(var_idxs, var_diffs))
    Alp.weighted_min_vertex_cover(m, distance)

    (Alp.get_option(m, :log_level) > 100) &&
        println("updated partition var selection => $(m.disc_vars)")
    return
end
