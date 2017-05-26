"""
    Set up a basic lifted mip model without any additional infromations.
    This model intends to be deep-copied at each iteration for a new MIP model.
    All non-linear term in the original model is lifted with a variable.

    This function is a prototype for mc-based bilinear model.
"""
function mcbi_bounding_mip(m::PODNonlinearModel; kwargs...)

    m.basic_model_mip = Model(solver=m.mip_solver)

    mcbi_post_basic_vars(m)
    mcbi_post_lifted_aff_constraints(m)
    mcbi_post_lifted_aff_obj(m)

    mcbi_post_mc_all(m)
end

function tightmccormick_one(m,xy,xz,y,xˡ,xᵘ,yˡ,yᵘ,z)
    @constraint(m, xy .>= xˡ*y + yˡ'*xz - xˡ*yˡ'*z)
    @constraint(m, xy .>= xᵘ*y + yᵘ'*xz - xᵘ*yᵘ'*z)
    @constraint(m, xy .<= xˡ*y + yᵘ'*xz - xˡ*yᵘ'*z)
    @constraint(m, xy .<= xᵘ*y + yˡ'*xz - xᵘ*yˡ'*z)
end

function tightmccormick_two(m,xy,xz,yz,xˡ,xᵘ,yˡ,yᵘ,zxy)
    @constraint(m, xy .>= xˡ'*yz + yˡ'*xz - xˡ'*zxy*yˡ)
    @constraint(m, xy .>= xᵘ'*yz + yᵘ'*xz - xᵘ'*zxy*yᵘ)
    @constraint(m, xy .<= xˡ'*yz + yᵘ'*xz - xˡ'*zxy*yᵘ)
    @constraint(m, xy .<= xᵘ'*yz + yˡ'*xz - xᵘ'*zxy*yˡ)
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
