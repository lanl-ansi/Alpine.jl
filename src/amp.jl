"""
    post_amp_convexification(m::PODNonlinearModel; kwargs...)

warpper function to convexify the problem for a bounding model. This function talks to nonlinear_info and convexification methods
to finish the last step required during the construction of bounding model.
"""
function post_amp_convexification(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    for i in 1:length(m.method_convexification)             # Additional user-defined convexificaition method
        eval(m.method_convexification[i])(m)
    end

    if !m.convex_disable_tmc
        post_amp_mccormick(m, use_discretization=discretization)    # handles all bi-linear and monomial convexificaitons
    end

    if !m.convex_disable_convhull
        # convex hull representation
    end

    convexification_exam(m)

    return
end

"""

"""
function post_amp_mccormick(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)

    # detect whether to use specific discretization information
    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    λ = Dict()
    λλ = Dict()
    λX = Dict()
    lb = Dict()
    ub = Dict()

    for bi in keys(m.nonlinear_info)
        if m.nonlinear_info[bi][:nonlinear_type] in [:monomial, :bilinear]
            @assert length(bi) == 2

            m.nonlinear_info[bi][:convexified] = true  # Bookeeping the examined terms
            idx_a = bi[1].args[2]
            idx_b = bi[2].args[2]
            idx_ab = m.nonlinear_info[bi][:lifted_var_ref].args[2]

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
                if m.nonlinear_info[bi][:nonlinear_type] == :monomial
                    mccormick_monomial(m.model_mip, Variable(m.model_mip, idx_ab), Variable(m.model_mip,idx_a), lb[idx_a][1], ub[idx_a][1])
                elseif m.nonlinear_info[bi][:nonlinear_type] == :bilinear
                    mccormick(m.model_mip, Variable(m.model_mip, idx_ab), Variable(m.model_mip, idx_a), Variable(m.model_mip, idx_b),
                        lb[idx_a][1], ub[idx_a][1], lb[idx_b][1], ub[idx_b][1])
                end
            else                                                    # Tighten McCormick
                # if m.nonlinear_info[bi][:monomial_satus]
                if m.nonlinear_info[bi][:nonlinear_type] == :monomial
                    λ = amp_post_λ(m.model_mip, λ, lb, ub, part_cnt_a, idx_a)
                    λX = amp_post_λX(m.model_mip, λX, part_cnt_a, idx_a, idx_b)
                    amp_post_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_a, idx_b)
                    amp_post_monomial_mc(m.model_mip, idx_ab, λ, λX, lb, ub, part_cnt_a, idx_a)
                # else
                elseif m.nonlinear_info[bi][:nonlinear_type] == :bilinear
                    # Partitioning on left
                    if (idx_a in m.var_discretization_mip) && !(idx_b in m.var_discretization_mip)
                        λ = amp_post_λ(m.model_mip, λ, lb, ub, part_cnt_a, idx_a)
                        λX = amp_post_λX(m.model_mip, λX, part_cnt_a, idx_a, idx_b)
                        λX[(idx_b,idx_a)] = [Variable(m.model_mip, idx_a)]
                        λλ = amp_post_λλ(m.model_mip, λλ, λ, idx_a, idx_b)
                        amp_post_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_a, idx_b)
                        amp_post_XX_mc(m.model_mip, idx_ab, λX, λλ, lb, ub, idx_a, idx_b)
                    end

                    # Partitioning of right
                    if !(idx_a in m.var_discretization_mip) && (idx_b in m.var_discretization_mip)
                        λ = amp_post_λ(m.model_mip, λ, lb, ub, part_cnt_b, idx_b)
                        λX = amp_post_λX(m.model_mip, λX, part_cnt_b, idx_b, idx_a)
                        λX[(idx_a,idx_b)] = [Variable(m.model_mip, idx_b)]
                        λλ = amp_post_λλ(m.model_mip, λλ, λ, idx_b, idx_a)
                        amp_post_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_b, idx_a)
                        amp_post_XX_mc(m.model_mip,idx_ab, λX, λλ, lb, ub, idx_b, idx_a)
                    end

                    # Partitioning on both variables
                    if (idx_a in m.var_discretization_mip) && (idx_b in m.var_discretization_mip)
                        λ = amp_post_λ(m.model_mip, λ, lb, ub, part_cnt_a, idx_a)
                        λ = amp_post_λ(m.model_mip, λ, lb, ub, part_cnt_b, idx_b)
                        λX = amp_post_λX(m.model_mip, λX, part_cnt_a, idx_a, idx_b)
                        λX = amp_post_λX(m.model_mip, λX, part_cnt_b, idx_b, idx_a)
                        λλ = amp_post_λλ(m.model_mip, λλ, part_cnt_a, part_cnt_b, idx_a, idx_b)
                        amp_post_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_a, idx_b)
                        amp_post_λxX_mc(m.model_mip, λX, λ, lb, ub, idx_b, idx_a)
                        amp_post_λxλ_mc(m.model_mip, λλ, λ, idx_a, idx_b)
                        amp_post_XX_mc(m.model_mip, idx_ab, λX, λλ, lb, ub, idx_a, idx_b)
                    end

                    # Error condition
                    if !(idx_a in m.var_discretization_mip) && !(idx_b in m.var_discretization_mip)
                        error("Error case. At least one term should show up in discretization choices.")
                    end
                end
            end
        end
    end
    return
end

function amp_post_λ(m::JuMP.Model, λ::Dict, lb::Dict, ub::Dict, dim::Int, idx::Int, relax=false)
    if !haskey(λ, idx)
        λ[idx] = @variable(m, [1:dim], Bin, basename=string("L",idx))
        @constraint(m, sum(λ[idx]) == 1) # The SOS-1 Constraints, not all MIP solver has SOS feature
        @constraint(m, Variable(m, idx) >= dot(lb[idx], λ[idx]))
        @constraint(m, Variable(m, idx) <= dot(ub[idx], λ[idx]))
    end
    return λ
end

function amp_post_monomial_mc(m::JuMP.Model, idx_aa::Int, λ::Dict, λX::Dict, LB::Dict, UB::Dict, dim::Int, idx_a::Int)
    @constraint(m, Variable(m, idx_aa) >= Variable(m, idx_a)^(2))
    @constraint(m, Variable(m, idx_aa) <= dot(λX[(idx_a,idx_a)],LB[idx_a]) + dot(λX[(idx_a,idx_a)],UB[idx_a]) - dot(λ[idx_a], *(diagm(LB[idx_a]), UB[idx_a])))
    return
end

function amp_post_λX(m::JuMP.Model, λX::Dict, dim::Int, idx_a::Int, idx_b::Int)
    if !haskey(λX, (idx_a,idx_b))
        λX[(idx_a,idx_b)] = @variable(m, [1:dim], basename=string("L",idx_a,"X",idx_b))
    end
    return λX
end

function amp_post_λλ(m::JuMP.Model, λλ::Dict, dim_a::Int, dim_b::Int, idx_a::Int, idx_b::Int)
    if !haskey(λλ, (idx_a, idx_b))
        λλ[(idx_a,idx_b)] = @variable(m, [1:dim_a, 1:dim_b], basename=string("L",idx_a,"L",idx_b))
        for i in 1:dim_a
    		for j in 1:dim_b
    			setlowerbound(λλ[(idx_a, idx_b)][i,j], 0);
    			setupperbound(λλ[(idx_a, idx_b)][i,j], 1);
    		end
    	end
        λλ[(idx_b,idx_a)] = λλ[(idx_a,idx_b)]'
    end
    return λλ
end

function amp_post_λλ(m::JuMP.Model, λλ::Dict, λ::Dict, idx_a::Int, idx_b::Int)
    if !haskey(λλ, (idx_a, idx_b))
        λλ[(idx_a,idx_b)] = λ[idx_a]
        λλ[(idx_a,idx_b)] = reshape(λλ[(idx_a,idx_b)], (length(λ[idx_a]), 1))
        λλ[(idx_b,idx_a)] = λλ[(idx_a,idx_b)]'
    end
    return λλ
end

function amp_post_XX_mc(m, ab, λX, λλ, LB, UB, a, b)
    @assert length(LB[a]) == length(UB[a])
    @assert length(LB[b]) == length(UB[b])
    dim_A = length(LB[a])
    dim_B = length(LB[b])
    @constraint(m, Variable(m, ab) .>= dot(λX[(a,b)],LB[a]) + dot(λX[(b,a)],LB[b]) - reshape(LB[a], (1, dim_A))*λλ[(a,b)]*reshape(LB[b], (dim_B, 1)))
	@constraint(m, Variable(m, ab) .>= dot(λX[(a,b)],UB[a]) + dot(λX[(b,a)],UB[b]) - reshape(UB[a], (1, dim_A))*λλ[(a,b)]*reshape(UB[b], (dim_B, 1)))
	@constraint(m, Variable(m, ab) .<= dot(λX[(a,b)],LB[a]) + dot(λX[(b,a)],UB[b]) - reshape(LB[a], (1, dim_A))*λλ[(a,b)]*reshape(UB[b], (dim_B, 1)))
	@constraint(m, Variable(m, ab) .<= dot(λX[(a,b)],UB[a]) + dot(λX[(b,a)],LB[b]) - reshape(UB[a], (1, dim_A))*λλ[(a,b)]*reshape(LB[b], (dim_B, 1)))
    return
end

function amp_post_λxX_mc(m::JuMP.Model, λX::Dict, λ::Dict, lb::Dict, ub::Dict, ind_λ::Int, ind_X::Int)

    # X_u and λ here are vectors, and X is one variable,
	# After the product, we have a polynomial to multiply X,
	# forming |λ| new variables stored in λX dictionary

	dim_λ = length(λ[ind_λ]) # This is how many new variables to be generated
	for i in 1:dim_λ
		lb_X = getlowerbound(Variable(m, ind_X))
		ub_X = getupperbound(Variable(m, ind_X))
		lb_λ = getlowerbound(λ[ind_λ][i])
		ub_λ = getupperbound(λ[ind_λ][i])
        @assert (lb_λ == 0.0) && (ub_λ == 1.0)
		mccormick(m, λX[(ind_λ,ind_X)][i], λ[ind_λ][i], Variable(m, ind_X), lb_λ, ub_λ, lb_X, ub_X)
	end
    return
end

function amp_post_λxλ_mc(m::JuMP.Model, λλ::Dict, λ::Dict, ind_A::Int, ind_B::Int)

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
			setlowerbound(λλ[(ind_A,ind_B)][i,j],0)
			setupperbound(λλ[(ind_A,ind_B)][i,j],1)
		end
	end

	for i in 1:dim_A
		for j in 1:dim_B
			mccormick_bin(m, λλ[(ind_A,ind_B)][i,j], λ[ind_A][i], λ[ind_B][j])
		end
	end

    return
end

function post_amp_vars(m::PODNonlinearModel; kwargs...)

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

function post_amp_lifted_constraints(m::PODNonlinearModel)

    for i in 1:m.num_constr_orig
        if m.constr_type_orig[i] == :(>=)
            @constraint(m.model_mip,
                sum(m.lifted_constr_aff_mip[i][:coefs][j]*Variable(m.model_mip, m.lifted_constr_aff_mip[i][:vars][j].args[2]) for j in 1:m.lifted_constr_aff_mip[i][:cnt]) >= m.lifted_constr_aff_mip[i][:rhs])
        elseif m.constr_type_orig[i] == :(<=)
            @constraint(m.model_mip,
                sum(m.lifted_constr_aff_mip[i][:coefs][j]*Variable(m.model_mip, m.lifted_constr_aff_mip[i][:vars][j].args[2]) for j in 1:m.lifted_constr_aff_mip[i][:cnt]) <= m.lifted_constr_aff_mip[i][:rhs])
        elseif m.constr_type_orig[i] == :(==)
            @constraint(m.model_mip,
                sum(m.lifted_constr_aff_mip[i][:coefs][j]*Variable(m.model_mip, m.lifted_constr_aff_mip[i][:vars][j].args[2]) for j in 1:m.lifted_constr_aff_mip[i][:cnt]) == m.lifted_constr_aff_mip[i][:rhs])
        end
    end

    return
end

function post_amp_lifted_objective(m::PODNonlinearModel)
    @objective(m.model_mip, m.sense_orig, sum(m.lifted_obj_aff_mip[:coefs][i]*Variable(m.model_mip, m.lifted_obj_aff_mip[:vars][i].args[2]) for i in 1:m.lifted_obj_aff_mip[:cnt]))
    return
end
