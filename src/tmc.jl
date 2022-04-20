function amp_post_mccormick(m::Optimizer; kwargs...)

    options = Dict(kwargs)

    # detect whether to use specific discretization information
    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization

    λ = Dict()
    λλ = Dict()
    λX = Dict()
    lb = Dict()
    ub = Dict()

    for bi in keys(m.nonconvex_terms)
        nl_type = m.nonconvex_terms[bi][:nonlinear_type]
        if ((!Alp.get_option(m, :monomial_convexhull))*(nl_type == :MONOMIAL) || (!Alp.get_option(m, :bilinear_convexhull))*(nl_type == :BILINEAR)) && (m.nonconvex_terms[bi][:convexified] == false)
            @assert length(bi) == 2
            m.nonconvex_terms[bi][:convexified] = true  # Bookkeeping the examined terms
            idx_a = bi[1].args[2]
            idx_b = bi[2].args[2]
            idx_ab = m.nonconvex_terms[bi][:lifted_var_ref].args[2]

            part_cnt_a = length(discretization[idx_a]) - 1
            part_cnt_b = length(discretization[idx_b]) - 1

            lb[idx_a] = discretization[idx_a][1:(end-1)]
            ub[idx_a] = discretization[idx_a][2:end]
            lb[idx_b] = discretization[idx_b][1:(end-1)]
            ub[idx_b] = discretization[idx_b][2:end]

            (lb[idx_a][1] == -Inf) && error("Infinite lower bound detected on VAR$idx_a, resolve it please")
            (lb[idx_b][1] == -Inf) && error("Infinite lower bound detected on VAR$idx_b, resolve it please")
            (ub[idx_a][end] == +Inf) && error("Infinite upper bound detected on VAR$idx_a, resolve it please")
            (ub[idx_b][end] == +Inf) && error("Infinite upper bound detected on VAR$idx_b, resolve it please")

            if (length(lb[idx_a]) == 1) && (length(lb[idx_b]) == 1)  # Basic McCormick
                if m.nonconvex_terms[bi][:nonlinear_type] == :MONOMIAL
                    Alp.mccormick_monomial(m.model_mip, _index_to_variable_ref(m.model_mip, idx_ab), _index_to_variable_ref(m.model_mip,idx_a), lb[idx_a][1], ub[idx_a][1])
                elseif m.nonconvex_terms[bi][:nonlinear_type] == :BILINEAR
                    Alp.mccormick(m.model_mip, _index_to_variable_ref(m.model_mip, idx_ab), _index_to_variable_ref(m.model_mip, idx_a), _index_to_variable_ref(m.model_mip, idx_b),
                        lb[idx_a][1], ub[idx_a][1], lb[idx_b][1], ub[idx_b][1])
                end
            else  # Tighten McCormick

                # if m.nonconvex_terms[bi][:MONOMIAL_satus]
                if m.nonconvex_terms[bi][:nonlinear_type] == :MONOMIAL
                    λ = Alp.amp_post_tmc_λ(m.model_mip, λ, lb, ub, part_cnt_a, idx_a)
                    λX = Alp.amp_post_tmc_λX(m.model_mip, λX, part_cnt_a, idx_a, idx_b)
                    Alp.amp_post_tmc_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_a, idx_b)
                    Alp.amp_post_tmc_monomial_mc(m.model_mip, idx_ab, λ, λX, lb, ub, part_cnt_a, idx_a)
                # else
                elseif m.nonconvex_terms[bi][:nonlinear_type] == :BILINEAR
                    # Partitioning on left
                    if (idx_a in m.disc_vars) && !(idx_b in m.disc_vars) && (part_cnt_b == 1)
                        λ = Alp.amp_post_tmc_λ(m.model_mip, λ, lb, ub, part_cnt_a, idx_a)
                        λX = Alp.amp_post_tmc_λX(m.model_mip, λX, part_cnt_a, idx_a, idx_b)
                        λX[(idx_b,idx_a)] = [_index_to_variable_ref(m.model_mip, idx_a)]
                        λλ = Alp.amp_post_tmc_λλ(λλ, λ, idx_a, idx_b)
                        Alp.amp_post_tmc_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_a, idx_b)
                        Alp.amp_post_tmc_XX_mc(m.model_mip, idx_ab, λX, λλ, lb, ub, idx_a, idx_b)
                    end

                    # Partitioning of right
                    if !(idx_a in m.disc_vars) && (idx_b in m.disc_vars) && (part_cnt_a == 1)
                        λ = Alp.amp_post_tmc_λ(m.model_mip, λ, lb, ub, part_cnt_b, idx_b)
                        λX = Alp.amp_post_tmc_λX(m.model_mip, λX, part_cnt_b, idx_b, idx_a)
                        λX[(idx_a,idx_b)] = [_index_to_variable_ref(m.model_mip, idx_b)]
                        λλ = Alp.amp_post_tmc_λλ(λλ, λ, idx_b, idx_a)
                        Alp.amp_post_tmc_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_b, idx_a)
                        Alp.amp_post_tmc_XX_mc(m.model_mip,idx_ab, λX, λλ, lb, ub, idx_b, idx_a)
                    end

                    # Partitioning on both variables
                    if (idx_a in m.disc_vars) && (idx_b in m.disc_vars)
                        λ = Alp.amp_post_tmc_λ(m.model_mip, λ, lb, ub, part_cnt_a, idx_a)
                        λ = Alp.amp_post_tmc_λ(m.model_mip, λ, lb, ub, part_cnt_b, idx_b)
                        λX = Alp.amp_post_tmc_λX(m.model_mip, λX, part_cnt_a, idx_a, idx_b)
                        λX = Alp.amp_post_tmc_λX(m.model_mip, λX, part_cnt_b, idx_b, idx_a)
                        λλ = Alp.amp_post_tmc_λλ(m.model_mip, λλ, part_cnt_a, part_cnt_b, idx_a, idx_b)
                        Alp.amp_post_tmc_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_a, idx_b)
                        Alp.amp_post_tmc_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_b, idx_a)
                        Alp.amp_post_tmc_λxλ_mc(m.model_mip, λλ, λ, idx_a, idx_b)
                        Alp.amp_post_tmc_XX_mc(m.model_mip, idx_ab, λX, λλ, lb, ub, idx_a, idx_b)
                    end

                    # Error condition
                    #if !(idx_a in m.disc_vars) && !(idx_b in m.disc_vars)
                    #    error("Error case. At least one term should show up in discretization choices.")
                    #end
                end
            end
        end
    end
    return
end

function amp_post_tmc_λ(m::JuMP.Model, λ::Dict, lb::Dict, ub::Dict, dim::Int, idx::Int, relax=false)
    if !haskey(λ, idx)
        λ[idx] = JuMP.@variable(m, [1:dim], Bin, base_name=string("L",idx))
        JuMP.@constraint(m, sum(λ[idx]) == 1) # SOS-1 Constraints
        JuMP.@constraint(m, _index_to_variable_ref(m, idx) >= dot(lb[idx], λ[idx]))
        JuMP.@constraint(m, _index_to_variable_ref(m, idx) <= dot(ub[idx], λ[idx]))
    end

    return λ
end

function amp_post_tmc_monomial_mc(m::JuMP.Model, idx_aa::Int, λ::Dict, λX::Dict, LB::Dict, UB::Dict, dim::Int, idx_a::Int)
    JuMP.@constraint(m, _index_to_variable_ref(m, idx_aa) >= _index_to_variable_ref(m, idx_a)^(2))
    JuMP.@constraint(m, _index_to_variable_ref(m, idx_aa) <= dot(λX[(idx_a,idx_a)],LB[idx_a]) + dot(λX[(idx_a,idx_a)],UB[idx_a]) - dot(λ[idx_a], *(Matrix(Diagonal(LB[idx_a])), UB[idx_a])))

    return
end

function amp_post_tmc_λX(m::JuMP.Model, λX::Dict, dim::Int, idx_a::Int, idx_b::Int)
    if !haskey(λX, (idx_a,idx_b))
        λX[(idx_a,idx_b)] = JuMP.@variable(m, [1:dim], base_name=string("L",idx_a,"X",idx_b))
    end

    return λX
end

function amp_post_tmc_λλ(m::JuMP.Model, λλ::Dict, dim_a::Int, dim_b::Int, idx_a::Int, idx_b::Int)
    if !haskey(λλ, (idx_a, idx_b))
        λλ[(idx_a,idx_b)] = JuMP.@variable(m, [1:dim_a, 1:dim_b], base_name=string("L",idx_a,"L",idx_b))
        for i in 1:dim_a
    		for j in 1:dim_b
    			JuMP.set_lower_bound(λλ[(idx_a, idx_b)][i,j], 0)
    			JuMP.set_upper_bound(λλ[(idx_a, idx_b)][i,j], 1)
    		end
    	end
        λλ[(idx_b,idx_a)] = λλ[(idx_a,idx_b)]'
    end

    return λλ
end

function amp_post_tmc_λλ(λλ::Dict, λ::Dict, idx_a::Int, idx_b::Int)
    if !haskey(λλ, (idx_a, idx_b))
        λλ[(idx_a,idx_b)] = λ[idx_a]
        λλ[(idx_a,idx_b)] = reshape(λλ[(idx_a,idx_b)], (length(λ[idx_a]), 1))
        λλ[(idx_b,idx_a)] = λλ[(idx_a,idx_b)]'
    end

    return λλ
end

function amp_post_tmc_XX_mc(m::JuMP.Model, ab, λX, λλ, LB, UB, a, b)
    @assert length(LB[a]) == length(UB[a])
    @assert length(LB[b]) == length(UB[b])
    dim_A = length(LB[a])
    dim_B = length(LB[b])
    JuMP.@constraint(m, _index_to_variable_ref(m, ab) .>= dot(λX[(a,b)],LB[a]) + dot(λX[(b,a)],LB[b]) - reshape(LB[a], (1, dim_A))*λλ[(a,b)]*reshape(LB[b], (dim_B, 1)))
	JuMP.@constraint(m, _index_to_variable_ref(m, ab) .>= dot(λX[(a,b)],UB[a]) + dot(λX[(b,a)],UB[b]) - reshape(UB[a], (1, dim_A))*λλ[(a,b)]*reshape(UB[b], (dim_B, 1)))
	JuMP.@constraint(m, _index_to_variable_ref(m, ab) .<= dot(λX[(a,b)],LB[a]) + dot(λX[(b,a)],UB[b]) - reshape(LB[a], (1, dim_A))*λλ[(a,b)]*reshape(UB[b], (dim_B, 1)))
	JuMP.@constraint(m, _index_to_variable_ref(m, ab) .<= dot(λX[(a,b)],UB[a]) + dot(λX[(b,a)],LB[b]) - reshape(UB[a], (1, dim_A))*λλ[(a,b)]*reshape(LB[b], (dim_B, 1)))

    return
end

_lower_bound(x) = JuMP.is_binary(x) ? 0.0 : JuMP.lower_bound(x)
_upper_bound(x) = JuMP.is_binary(x) ? 1.0 : JuMP.upper_bound(x)

function amp_post_tmc_λxX_mc(m::JuMP.Model, λX::Dict, λ::Dict, lb::Dict, ub::Dict, ind_λ::Int, ind_X::Int)

    # X_u and λ here are vectors, and X is one variable,
	# After the product, we have a polynomial to multiply X,
	# forming |λ| new variables stored in λX dictionary

	dim_λ = length(λ[ind_λ]) # This is how many new variables to be generated
	for i in 1:dim_λ
        v = _index_to_variable_ref(m, ind_X)
		lb_X = JuMP.lower_bound(v)
		ub_X = JuMP.upper_bound(v)
        lb_λ = _lower_bound(λ[ind_λ][i])
		ub_λ = _upper_bound(λ[ind_λ][i])
        @assert (lb_λ == 0.0) && (ub_λ == 1.0)
		Alp.mccormick(m, λX[(ind_λ,ind_X)][i], λ[ind_λ][i], v, lb_λ, ub_λ, lb_X, ub_X)
	end

    return
end

function amp_post_tmc_λxλ_mc(m::JuMP.Model, λλ::Dict, λ::Dict, ind_A::Int, ind_B::Int)

	#=
    A 3 x 2 example, some new variables are formed
	| y11y21 y11y22 |   | y11 |   | y21 |T
	| y12y21 y12y22 | = | y12 | x |     |
	| y13y21 y13y22 |   | y13 |	  |	y22 |
	=#

	# Prepare the new variable, these two counts track the total count of new variables
	dim_A = length(λ[ind_A])
	dim_B = length(λ[ind_B])
	for i in 1:dim_A
		for j in 1:dim_B
			JuMP.set_lower_bound(λλ[(ind_A,ind_B)][i,j],0)
			JuMP.set_upper_bound(λλ[(ind_A,ind_B)][i,j],1)
		end
	end

	for i in 1:dim_A
		for j in 1:dim_B
			Alp.mccormick_bin(m, λλ[(ind_A,ind_B)][i,j], λ[ind_A][i], λ[ind_B][j])
		end
	end

    return
end

"""
    mccormick(m::JuMP.Model, xy, x, y, x_l, x_u, y_l, y_u)

Generic function to add a McCormick convex envelop, where `xy=x*y` and `x_l, x_u, y_l, y_u` are variable bounds.
"""
function mccormick(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, y::JuMP.VariableRef, lb_x, ub_x, lb_y, ub_y)

    JuMP.@constraint(m, xy >= lb_x*y + lb_y*x - lb_x*lb_y)
    JuMP.@constraint(m, xy >= ub_x*y + ub_y*x - ub_x*ub_y)
    JuMP.@constraint(m, xy <= lb_x*y + ub_y*x - lb_x*ub_y)
    JuMP.@constraint(m, xy <= ub_x*y + lb_y*x - ub_x*lb_y)

    return
end

function mccormick_binlin(m::JuMP.Model, binlin::JuMP.VariableRef, bin::JuMP.VariableRef, lin::JuMP.VariableRef, lb, ub)

    # TODO think about how to address this issue
    warnuser = false
    if ub == Inf
        ub = 1e6
        warnuser = true
    end

    if lb == -Inf
        lb = -1e6
        warnuser = true
    end

    if lb >= 0
        JuMP.@constraint(m, binlin <= ub*bin)
        JuMP.@constraint(m, binlin <= lin)
        JuMP.@constraint(m, binlin >= lin - (1-bin)*ub)
    else
        JuMP.@constraint(m, binlin <= ub)
        JuMP.@constraint(m, binlin >= lb)
        JuMP.@constraint(m, binlin <= bin*ub)
        JuMP.@constraint(m, binlin >= bin*lb)
        JuMP.@constraint(m, binlin <= lin - (1-bin)*lb)
        JuMP.@constraint(m, binlin >= lin - (1-bin)*ub)
    end

    # Second position to handle inf bounds
    warnuser && @warn "BINLIN term exception using -1e6/1e6 as lb/ub"

    return
end

function mccormick_bin(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, y::JuMP.VariableRef)
    JuMP.@constraint(m, xy <= x)
    JuMP.@constraint(m, xy <= y)
    JuMP.@constraint(m, xy >= x+y-1)
    return
end

# [TODO] Unused functions but will be used for later
function mccormick_monomial(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, lb_x, ub_x)
    JuMP.@constraint(m, xy >= x^2)
    JuMP.@constraint(m, xy <= (lb_x+ub_x)*x - (lb_x*ub_x))

    return
end

# Fortet linearization
# Reference: https://doi.org/10.1007/s10288-006-0015-3
function binprod_relax(m::JuMP.Model, z::JuMP.VariableRef, x::Vector)

    for i in x
        JuMP.@constraint(m, z <= i)
    end
    JuMP.@constraint(m, z >= sum(x) - (length(x)-1))

    return
end

#=
function tightmccormick_monomial(m,x_p,x,xz,lb_x,ub_x,z,p,lazy,quad) # if p=2, tightened_lazycuts = tightmccormick_quad
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
        JuMP.@constraint(m, x_p >= x^2)
    else
        x0_vec = sort(union(lb_x, ub_x))
        for x0 in x0_vec
            JuMP.@constraint(m, x_p >= (1-p)*(x0)^p + p*(x0)^(p-1)*x)
        end
    end

    A = ((ub_x).^p-(lb_x).^p)./(ub_x-lb_x)
    JuMP.@constraint(m, x_p .<= A'*xz - (A.*lb_x)'*z + ((lb_x).^p)'*z)

    return
end
=#
