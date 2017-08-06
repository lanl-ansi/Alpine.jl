"""

    create_bounding_mip(m::PODNonlinearModel; use_discretization::Dict)

Set up a JuMP MILP bounding model base on variable domain partitioning information stored in `use_discretization`.
By default, if `use_discretization is` not provided, it will use `m.discretizations` store in the POD model.
The basic idea of this MILP bounding model is to use Tighten McCormick to convexify the original Non-convex region.
Among all presented partitionings, the bounding model will choose one specific partition as the lower bound solution.
The more partitions there are, the better or finer bounding model relax the original MINLP while the more
efforts required to solve this MILP is required.

This function is implemented in the following manner:

    * [`amp_post_vars`](@ref): post original and lifted variables
    * [`amp_post_lifted_constraints`](@ref): post original and lifted constraints
    * [`amp_post_lifted_obj`](@ref): post original or lifted objective function
    * [`amp_post_tmc_mccormick`](@ref): post Tighen McCormick variables and constraints base on `discretization` information

More specifically, the Tightening McCormick used here can be genealized in the following mathematcial formulation. Consider a nonlinear term
```math
\\begin{subequations}
\\begin{align}
   &\\widehat{x_{ij}} \\geq (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\widehat{x_{ij}} \\geq (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\widehat{x_{ij}} \\leq (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\widehat{x_{ij}} \\leq (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i) x_j + (\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) x_i - (\\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i)(\\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) \\\\
   & \\mathbf{x}_i^u\\cdot\\hat{\\mathbf{y}}_i) \\geq x_{i} \\geq \\mathbf{x}_i^l\\cdot\\hat{\\mathbf{y}}_i) \\\\
   & \\mathbf{x}_j^u\\cdot\\hat{\\mathbf{y}}_j) \\geq x_{j} \\geq \\mathbf{x}_j^l\\cdot\\hat{\\mathbf{y}}_j) \\\\
   &\\sum \\hat{\\mathbf{y}_i} = 1, \\ \\ \\sum \\hat{\\mathbf{y}_j}_k = 1 \\\\
   &\\hat{\\mathbf{y}}_i \\in \\{0,1\\}, \\hat{\\mathbf{y}}_j \\in \\{0,1\\}
\\end{align}
\\end{subequations}
```

"""
function create_bounding_mip(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP Model
    start_build = time()
    # ------- Model Construction ------ #
    amp_post_vars(m)                                                # Post original and lifted variables
    amp_post_lifted_constraints(m)                                  # Post lifted constraints
    amp_post_lifted_objective(m)                                    # Post objective
    amp_post_convexification(m, use_discretization=discretization)  # Convexify problem
    # --------------------------------- #
    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

    return
end

"""
    amp_post_convexification(m::PODNonlinearModel; kwargs...)

warpper function to convexify the problem for a bounding model. This function talks to nonlinear_terms and convexification methods
to finish the last step required during the construction of bounding model.
"""
function amp_post_convexification(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    for i in 1:length(m.method_convexification)             # Additional user-defined convexification method
        eval(m.method_convexification[i])(m)
    end

    if !m.bilinear_mccormick
        amp_post_tmc_mccormick(m, use_discretization=discretization)    # handles all bi-linear and monomial convexificaitons
    end

    if !m.bilinear_convexhull
        # convex hull representation
        # amp_post_convhull(m, use_discretization=discretization)
    end

    convexification_exam(m)

    return
end

function amp_post_vars(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    if haskey(options, :use_discretization)
        l_var = [options[:use_discretization][i][1]   for i in 1:(m.num_var_orig+m.num_var_lifted_mip)]
        u_var = [options[:use_discretization][i][end] for i in 1:(m.num_var_orig+m.num_var_lifted_mip)]
    else
        l_var = m.l_var_tight
        u_var = m.u_var_tight
    end

    @variable(m.model_mip, x[i=1:(m.num_var_orig+m.num_var_lifted_mip)])
    for i in 1:(m.num_var_orig+m.num_var_lifted_mip)
        (i <= m.num_var_orig) && setcategory(x[i], m.var_type_orig[i])
        (l_var[i] > -Inf) && (setlowerbound(x[i], l_var[i]))   # Changed to tight bound, if no bound tightening is performed, will be just .l_var_orig
        (u_var[i] < Inf) && (setupperbound(x[i], u_var[i]))   # Changed to tight bound, if no bound tightening is performed, will be just .u_var_orig
    end

    return
end

function amp_post_lifted_constraints(m::PODNonlinearModel)

    for i in 1:m.num_constr_orig
        if m.structural_constr[i] == :affine
            amp_post_affine_constraint(m.model_mip, m.bounding_constr_mip[i])
        elseif m.structural_constr[i] == :convex
            amp_post_convex_constraint(m.model_mip, m.bounding_constr_mip[i])
        else
            error("Unknown structural_constr type $(m.structural_constr[i])")
        end
    end

    return
end

function amp_post_affine_constraint(model_mip::JuMP.Model, affine::Dict)

    if affine[:sense] == :(>=)
        @constraint(model_mip,
            sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) >= affine[:rhs])
    elseif affine[:sense] == :(<=)
        @constraint(model_mip,
            sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) <= affine[:rhs])
    elseif affine[:sense] == :(==)
        @constraint(model_mip,
            sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) == affine[:rhs])
    end

    return
end

function amp_post_convex_constraint(model_mip::JuMP.Model, convex::Dict)

    if convex[:sense] == :(>=)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) >= convex[:rhs])
    elseif convex[:sense] == :(<=)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) <= convex[:rhs])
    elseif convex[:sense] == :(==)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) == convex[:rhs])
    end

    return
end

function amp_post_lifted_objective(m::PODNonlinearModel)

    if m.structural_obj == :affine
        @objective(m.model_mip, m.sense_orig, m.bounding_obj_mip[:rhs] + sum(m.bounding_obj_mip[:coefs][i]*Variable(m.model_mip, m.bounding_obj_mip[:vars][i].args[2]) for i in 1:m.bounding_obj_mip[:cnt]))
    elseif m.structural_obj == :convex
        @objective(m.model_mip, m.sense_orig, m.bounding_obj_mip[:rhs] + sum(m.bounding_obj_mip[:coefs][i]*Variable(m.model_mip, m.bounding_obj_mip[:vars][i].args[2])^2 for i in 1:m.bounding_obj_mip[:cnt]))
    else
        error("Unknown structural obj type $(m.structural_obj)")
    end

    return
end
