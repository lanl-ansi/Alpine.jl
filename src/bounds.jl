@doc """
    Initialize the discretization of each variables by setting up the bounds.
""" ->
function initialize_discretization(m::PODNonlinearModel;kwargs...)
    for var in 1:m.num_var_orig
        # TODO :: add a scheme for equal bounds
        lb = m.l_var_orig[var]
        ub = m.u_var_orig[var]
        m.discretization[var] = [lb, ub]
    end

    return
end

# function initialize_discretization(m::PODNonlinearModel; kwargs...)
#     @printf "running algorithm for choosing the variables to discretize\n"
#     pick_vars_discretization(m)
#     for var in 1:m.num_var_orig
#         lb = m.l_var_orig[var]
#         ub = m.u_var_orig[var]
#         if var in m.var_discretization_mip
#             # m.discretization[var] = [lb, m.sol_incumb_ub[var], ub]  # Alternative way to construct
#             point = m.best_sol[var]
#             radius = (ub-lb)/m.discretization_ratio
#             local_lb = max(lb, point-radius)
#             local_ub = min(ub, point+radius)
#             m.discretization[var] = unique([lb, local_lb, local_ub, ub])
#         else
#             m.discretization[var] = [lb, ub]
#         end
#     end
# end
