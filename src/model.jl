"""
    Set up a basic lifted mip model without any additional infromations.
    This model intends to be deep-copied at each iteration for a new MIP model. 
    All non-linear term in the original model is lifted with a variable.
    The performance 
"""
function setup_basic_bounding_mip(m::PODNonlinearModel; kwargs...)

    # Setup the MIP model and hook the solver
    m.basic_model_mip = Model(solver=m.mip_solver)

    # Setup original variables (following original setup, warmstart omitted not utilized)
    @variable(m.basic_model_mip, m.l_var_orig[i]<=x[i=1:m.num_var_orig]<=m.u_var_orig[i], category=m.var_type_orig[i])

    # Setup lifted variables (no specific setup)
    @variable(m.basic_model_mip, x_lift[(m.num_var_orig+1):(m.num_var_orig+m.lifted_x_cont)])

    # Setup Constraints Function (lifted if necessary)
    for i in 1:m.num_constr_orig
        ex = deepcopy(m.lifted_constr_expr_mip[i])  # Make a copy Keep the fiedility of the lifted expression
        expr_dereferencing(ex.args[2], m.basic_model_mip) # Dereferencing is in-place and model specific
        JuMP.addNLconstraint(m.basic_model_mip, ex) # Using addNLconstraint since it is :Expr friendly
    end

    # Setup Objective Function (lifted if necessary)
    ex = deepcopy(m.lifted_obj_expr_mip)
    expr_dereferencing(ex, m.basic_model_mip)
    JuMP.setNLobjective(m.basic_model_mip, m.sense_orig, ex) # Using addNLobjective since it is :Expr friendly

    return m
end

