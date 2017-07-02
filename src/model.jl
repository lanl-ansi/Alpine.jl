"""

    create_bounding_mip(m::PODNonlinearModel; use_discretization::Dict)

Set up a JuMP MILP bounding model base on variable domain partitioning information stored in `use_discretization`.
By default, if `use_discretization is` not provided, it will use `m.discretizations` store in the POD model.
The basic idea of this MILP bounding model is to use Tighten McCormick to convexify the original Non-convex region.
Among all presented partitionings, the bounding model will choose one specific partition as the lower bound solution.
The more partitions there are, the better or finer bounding model relax the original MINLP while the more
efforts required to solve this MILP is required.

This function is implemented in the following manner:

    * [`post_amp_vars`](@ref): post original and lifted variables
    * [`post_amp_lifted_constraints`](@ref): post original and lifted constraints
    * [`post_amp_lifted_obj`](@ref): post original or lifted objective function
    * [`post_amp_mccormick`](@ref): post Tighen McCormick variables and constraints base on `discretization` information

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
    post_amp_vars(m)                                            # Post original and lifted variables
    post_amp_lifted_constraints(m)                              # Post lifted constraints
    post_amp_lifted_objective(m)                                # Post objective
    post_amp_mccormick(m, use_discretization=discretization)    # handles all bi-linear and monomial convexificaitons
    # any additional built-in convexification method goes here
    for i in 1:length(m.expression_convexification)             # Additional user-defined convexificaition method
        eval(m.expression_convexification[i])(m)
    end
    # --------------------------------- #
    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

    return
end

"""
    pick_vars_discretization(m::PODNonlinearModel)

This function helps pick the variables for discretization. The method chosen depends on user-inputs.
In case when `indices::Int` is provided, the method is chosen as built-in method. Currently,
there exist two built-in method:

    * `max-cover(m.discretization_var_pick_algo=0, default)`: pick all variables involved in the non-linear term for discretization
    * `min-vertex-cover(m.discretization_var_pick_algo=1)`: pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified

For advance usage, `m.discretization_var_pick_algo` allows `::Function` inputs. User is required to perform flexible methods in choosing the non-linear variable.
For more information, read more details at [Hacking Solver](@ref).

"""
function pick_vars_discretization(m::PODNonlinearModel)

    if isa(m.discretization_var_pick_algo, Function)
        (m.log_level > 0) && println("method for picking discretization variable provided... overiding the parameter input...")
        eval(m.discretization_var_pick_algo)(m)
        # TODO: Perform generic validation for user-defined method
    elseif isa(m.discretization_var_pick_algo, Int)
        # Figure out which are the variables that needs to be partitioned
        if m.discretization_var_pick_algo == 0
            max_cover(m)
        elseif m.discretization_var_pick_algo == 1
            min_vertex_cover(m)
        else
            error("Unsupported default method for picking variables for discretization")
        end
    else
        error("Input for parameter :discretization_var_pick_algo is illegal. Should be either a Int for default methods indexes or functional inputs.")
    end

    return
end


function mccormick(m,xy,x,y,xˡ,xᵘ,yˡ,yᵘ)

    @constraint(m, xy >= xˡ*y + yˡ*x - xˡ*yˡ)
    @constraint(m, xy >= xᵘ*y + yᵘ*x - xᵘ*yᵘ)
    @constraint(m, xy <= xˡ*y + yᵘ*x - xˡ*yᵘ)
    @constraint(m, xy <= xᵘ*y + yˡ*x - xᵘ*yˡ)

    return
end

function mccormick_bin(m,xy,x,y)

    @constraint(m, xy <= x)
    @constraint(m, xy <= y)
    @constraint(m, xy >= x+y-1)

    return
end

function mccormick_monomial(m,xy,x,xˡ,xᵘ)
    @constraint(m, xy >= x^2)
    @constraint(m, xy <= (xˡ+xᵘ)*x - (xˡ*xᵘ))

    return
end

function tightmccormick_monomial(m,x_p,x,xz,xˡ,xᵘ,z,p,lazy,quad) # if p=2, tightened_lazycuts = tightmccormick_quad
    if lazy == 1
        function GetLazyCuts_quad(cb)
            TOL1 = 1e-6
            if (getvalue(x)^p > (getvalue(x_p) + TOL1))
                a = p*getvalue(x)^(p-1)
                b = (1-p)*getvalue(x)^p
                @lazyconstraint(cb, a*x + b <= x_p)
            end
        end
        addlazycallback(m, GetLazyCuts_quad)
    elseif p == 2 && quad == 1
        @constraint(m, x_p >= x^2)
    else
        x0_vec = sort(union(xˡ, xᵘ))
        for x0 in x0_vec
            @constraint(m, x_p >= (1-p)*(x0)^p + p*(x0)^(p-1)*x)
        end
    end

    A = ((xᵘ).^p-(xˡ).^p)./(xᵘ-xˡ)
    @constraint(m, x_p .<= A'*xz - (A.*xˡ)'*z + ((xˡ).^p)'*z)

    return
end

function post_obj_bounds(m::PODNonlinearModel, bound::Float64; kwargs...)
    if m.sense_orig == :Max
        @constraint(m.model_mip,
            sum(m.lifted_obj_aff_mip[:coefs][j]*Variable(m.model_mip, m.lifted_obj_aff_mip[:vars][j].args[2]) for j in 1:m.lifted_obj_aff_mip[:cnt]) >= bound)
    elseif m.sense_orig == :Min
        @constraint(m.model_mip,
            sum(m.lifted_obj_aff_mip[:coefs][j]*Variable(m.model_mip, m.lifted_obj_aff_mip[:vars][j].args[2]) for j in 1:m.lifted_obj_aff_mip[:cnt]) <= bound)
    end

    return
end

"""

    min_vertex_cover(m:PODNonlinearModel)

A built-in method for selecting variables for discretization.

"""
function min_vertex_cover(m::PODNonlinearModel)

    # Collect the information for arcs and nodes
    nodes = Set()
    arcs = Set()
    for pair in keys(m.nonlinear_info)
        arc = []
        for i in pair
            @assert isa(i.args[2], Int)
            push!(nodes, i.args[2])
            push!(arc, i.args[2])
        end
        push!(arcs, arc)
    end
    nodes = collect(nodes)
    arcs = collect(arcs)

    # Set up minimum vertex cover problem
    minvertex = Model(solver=m.mip_solver)
    @variable(minvertex, x[nodes], Bin)
    for arc in arcs
        @constraint(minvertex, x[arc[1]] + x[arc[2]] >= 1)
    end
    @objective(minvertex, Min, sum(x))
    status = solve(minvertex, suppress_warnings=true)

    xVal = getvalue(x)

    # Getting required information
    m.num_var_discretization_mip = Int(sum(xVal))
    m.var_discretization_mip = [i for i in nodes if xVal[i] > 1e-5]

end

"""

    max_cover(m:PODNonlinearModel)

A built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.

"""
function max_cover(m::PODNonlinearModel; kwargs...)

    nodes = Set()
    for pair in keys(m.nonlinear_info)
        for i in pair
            @assert isa(i.args[2], Int)
            push!(nodes, i.args[2])
        end
    end
    nodes = collect(nodes)
    m.num_var_discretization_mip = length(nodes)
    m.var_discretization_mip = nodes

    return
end
