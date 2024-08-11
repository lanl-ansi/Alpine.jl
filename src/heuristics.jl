"""
    Ranking of variables involved in nonlinear terms for piecewise adaptive partitioning:
    Ranked based on the absolute gap between each variable's solution from the lower-bounding MIP and the best feasible solution to the MINLP.
    Currently doesn't support recursive convexification
"""
function update_disc_cont_var(m::Optimizer)
    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    var_idxs = copy(m.candidate_disc_vars)
    var_diffs = Vector{Float64}(
        undef,
        m.num_var_orig + length(keys(m.linear_terms)) + length(keys(m.nonconvex_terms)),
    )

    for i in 1:m.num_var_orig       # Original Variables
        var_diffs[i] = abs(m.best_sol[i] - m.best_bound_sol[i])
    end

    !isempty(m.linear_terms) && (var_diffs = _eval_var_diffs!(m.linear_terms, var_diffs))
    !isempty(m.nonconvex_terms) &&
        (var_diffs = _eval_var_diffs!(m.nonconvex_terms, var_diffs))

    distance = Dict(zip(var_idxs, var_diffs))
    weighted_min_vertex_cover(m, distance)

    (get_option(m, :log_level) > 100) &&
        println("updated partition var selection => $(m.disc_vars)")
    return
end

# Helper function for sequential evaluation to avoid dependency issue
function _eval_var_diffs!(terms::Dict{Any,Any}, var_diffs::Vector{Float64})
    for i in 1:length(keys(terms))
        for j in keys(terms)
            if terms[j][:id] == i
                var_diffs[terms[j][:lifted_var_ref].args[2]] =
                    terms[j][:evaluator](terms[j], var_diffs)
            end
        end
    end

    return var_diffs
end
