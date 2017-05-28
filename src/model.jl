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
    return m
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
    pick_discretize_vars(m)
    for var in 1:m.num_var_orig
        if var in m.discrete_x
            m.discretization[var] = [m.l_var_orig[var], m.l_var_orig[var]+(m.u_var_orig[var]-m.l_var_orig[var])/2, m.u_var_orig[var]]
        else
            m.discretization[var] = [m.l_var_orig[var], m.u_var_orig[var]]
        end
    end
end
