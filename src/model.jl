"""
    Set up a basic lifted mip model without any additional infromations.
    This model intends to be deep-copied at each iteration for a new MIP model.
    All non-linear term in the original model is lifted with a variable.

    This function is a prototype for mc-based bilinear model.
"""
function lower_bounding_mip(m::PODNonlinearModel; kwargs...)

    m.basic_model_mip = Model(solver=m.mip_solver) # Construct model

    start_lb_build = time()
    mcbi_post_basic_vars(m)             # Post original and lifted variables
    mcbi_post_lifted_aff_constraints(m) # Post lifted constraints
    mcbi_post_lifted_aff_obj(m)         # Post objective
    mcbi_post_mc(m)                     # Post all mccormick problem
    cputime_lb_build = time() - start_lb_build
    m.logs[:total_time] += cputime_lb_build

end

function pick_discretize_vars(m::PODNonlinearModel)
    # Figure out which are the variables that needs to be partitioned
    if m.discrete_vars_choice == 0
        max_cover(m)
    elseif m.discrete_vars_choice == 1
        min_vertex_cover(m)
    else
        error("Unsupported method for picking variables for discretization")
    end
end

function mccormick(m,xy,x,y,xˡ,xᵘ,yˡ,yᵘ)
    @constraint(m, xy >= xˡ*y + yˡ*x - xˡ*yˡ)
    @constraint(m, xy >= xᵘ*y + yᵘ*x - xᵘ*yᵘ)
    @constraint(m, xy <= xˡ*y + yᵘ*x - xˡ*yᵘ)
    @constraint(m, xy <= xᵘ*y + yˡ*x - xᵘ*yˡ)
end

function mccormick_bin(m,xy,x,y)
    @constraint(m, xy <= x)
    @constraint(m, xy <= y)
    @constraint(m, xy >= x+y-1)
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
    return m
end

function initialize_discretization(m::PODNonlinearModel; kwargs...)
    options = Dict(kwargs)
    pick_discretize_vars(m)
    for var in 1:m.num_var_orig
        lb = m.l_var_orig[var]
        ub = m.u_var_orig[var]
        if var in m.discrete_x
            # m.discretization[var] = [lb, m.sol_incumb_ub[var], ub]  # Alternative way to construct
            point = m.sol_incumb_ub[var]
            radius = (ub-lb)/m.discrete_ratio
            local_lb = max(lb, point-radius)
            local_ub = min(ub, point+radius)
            m.discretization[var] = unique([lb, local_lb, local_ub, ub])
        else
            m.discretization[var] = [lb, ub]
        end
    end
end

"""
    Notice there are other ways to construct partitions
"""
function add_discretization(m::PODNonlinearModel; kwargs...)

    #=
    Consider original partition [0, 3, 7, 9], where LB solution is 4.
    Use ^ as the new partition, | as the original partition

    A case when discretize ratio = 8
    | -------- | - ^ -- * -- ^ ---- | -------- |
    0          3  3.5   4   4.5     7          9

    A special case when discretize ratio = 4
    | -------- | --- * --- ^ ---- | -------- |
    0          3     4     5      7          9
    =#

    for i in 1:m.num_var_orig
        point = m.sol_incumb_lb[i]
        if i in m.discrete_x  # Only construct when discretized
            # @show "Before ", i, m.discretization[i]
            for j in 1:length(m.discretization[i])
                if point >= m.discretization[i][j] && point <= m.discretization[i][j+1]  # Locating the right location
                    @assert j < length(m.discretization[i])
                    lb_local = m.discretization[i][j]
                    ub_local = m.discretization[i][j+1]
                    distance = ub_local - lb_local
                    radius = distance / m.discrete_ratio
                    lb_new = max(point - radius/2, lb_local)
                    ub_new = min(point + radius/2, ub_local)
                    # @show j, point, lb_new, ub_new, lb_local, ub_local
                    if ub_new < ub_local  # Insert new UB-based partition
                        insert!(m.discretization[i], j+1, ub_new)
                    end
                    if lb_new > lb_local  # Insert new LB-based partition
                        insert!(m.discretization[i], j+1, lb_new)
                    end
                    break
                end
            end
            # @show "After ", i, m.discretization[i], "\n"
        end
    end

end


"""
    Use lower bound solution active interval to tight upper bound variable bounds
"""
function tight_ub_bounds(m::PODNonlinearModel; kwargs...)

    l_var = copy(m.l_var_orig)
    u_var = copy(m.u_var_orig)
    for i in 1:m.num_var_orig
        if i in m.discrete_x
            point = m.sol_incumb_lb[i]
            for j in 1:length(m.discretization[i])
                if point >= m.discretization[i][j] && point <= m.discretization[i][j+1]
                    @assert j < length(m.discretization[i])
                    l_var[i] = m.discretization[i][j]
                    u_var[i] = m.discretization[i][j+1]
                end
            end
        end
    end
    return l_var, u_var
end


"""
    Perform a minimum vertex cover for selecting variables for discretization
"""
function min_vertex_cover(m::PODNonlinearModel; kwargs...)

    # Collect the information for arcs and nodes
    nodes = Set()
    arcs = Set()
    for pair in keys(m.dict_nonlinear_info)
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
    m.discrete_x_cnt = Int(sum(xVal))
    m.discrete_x = [i for i in nodes if xVal[i]>0]

end

"""
    Collect all envolved variables as a max cover for generate discretization
"""
function max_cover(m::PODNonlinearModel; kwargs...)
    nodes = Set()
    for pair in keys(m.dict_nonlinear_info)
        for i in pair
            @assert isa(i.args[2], Int)
            push!(nodes, i.args[2])
        end
    end
    nodes = collect(nodes)
    m.discrete_x_cnt = length(nodes)
    m.discrete_x = nodes
end
