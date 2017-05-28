"""
    Set up a basic lifted mip model without any additional infromations.
    This model intends to be deep-copied at each iteration for a new MIP model.
    All non-linear term in the original model is lifted with a variable.

    This function is a prototype for mc-based bilinear model.
"""
function mcbi_bounding_mip(m::PODNonlinearModel; kwargs...)

    m.basic_model_mip = Model(solver=m.mip_solver) # Construct model

    mcbi_post_basic_vars(m)             # Post original and lifted variables
    mcbi_post_lifted_aff_constraints(m) # Post lifted constraints
    mcbi_post_lifted_aff_obj(m)         # Post objective
    mcbi_post_mc(m)                     # Post all mccormick problem

end

"""
    Compact function for all partioning variables. High-level wrapper for constructing all mccormick convexifications.
"""
function mcbi_post_mc(m::PODNonlinearModel; kwargs...)

    λ = Dict()
    λλ = Dict()
    λX = Dict()
    lb = Dict()
    ub = Dict()

    for bi in keys(m.dict_nonlinear_info)

        # Consider a general bilinear framework c = a * b
        # part -> short for partitions

        @assert length(bi) == 2
        idx_a = bi[1].args[2]
        idx_b = bi[2].args[2]
        idx_ab = m.dict_nonlinear_info[bi][:lifted_var_ref].args[2]

        part_cnt_a = length(m.discretization[idx_a]) - 1
        part_cnt_b = length(m.discrsetization[idx_b]) - 1

        lb[idx_a] = m.discretization[idx_a][1:(end-1)] # Hope Julia can link the list efficiently
        ub[idx_a] = m.discretization[idx_a][2:end]
        lb[idx_b] = m.discretization[idx_b][1:(end-1)]
        ub[idx_b] = m.discretization[idx_b][2:end]

        # Partitioning on left
        if (idx_a in m.discrete_x) && !(idx_b in m.discrete_x)
            λ = mcbi_post_λ(m.basic_model_mip, λ, part_cnt_a, idx_a)
            λX = mcbi_post_λX(m.basic_model_mip, λX, part_cnt_a, idx_a, idx_b)
            λX[(idx_b,idx_a)] = [Variable(m.basic_model_mip, idx_b)]
            λλ = mcbi_post_λλ(m.basic_model_mip, λλ, λ, idx_a, idx_b)
            mcbi_post_λxX_mc(m.basic_model_mip, λX, λ, idx_a, idx_b)
            mcbi_post_XX_mc(m.basic_model_mip, idx_ab, λX, λλ, lb, ub, idx_a, idx_b)
        end

        # Partitioning of right
        if !(idx_a in m.discrete_x) && (idx_b in m.discrete_x)
            λ = mcbi_post_λ(m.basic_model_mip, λ, part_cnt_b, idx_b)
            λX = mcbi_post_λX(m.basic_model_mip, λX, part_cnt_b, idx_b, idx_a)
            λX[(idx_a,idx_b)] = [Variable(m.basic_model_mip, idx_a)]
            λλ = mcbi_post_λλ(m.basic_model_mip, λλ, λ, idx_b, idx_a)
            mcbi_post_λxX_mc(m.basic_model_mip, λX, λ, idx_b, idx_a)
            mcbi_post_XX_mc(m.basic_model_mip,idx_ab, λX, λλ, lb, ub, idx_b, idx_a)
        end

        # Partitioning on both variables
        if (idx_a in m.discrete_x) && (idx_b in m.discrete_x)
            λ = mcbi_post_λ(m.basic_model_mip, λ, part_cnt_a, idx_a)
            λ = mcbi_post_λ(m.basic_model_mip, λ, part_cnt_b, idx_b)
            λX = mcbi_post_λX(m.basic_model_mip, λX, part_cnt_a, idx_a, idx_b)
            λX = mcbi_post_λX(m.basic_model_mip, λX, part_cnt_b, idx_b, idx_a)
            λλ = mcbi_post_λλ(m.basic_model_mip, λλ, part_cnt_a, part_cnt_b, idx_a, idx_b)
            mcbi_post_λxX_mc(m.basic_model_mip, λX, λ, idx_a, idx_b)
            mcbi_post_λxX_mc(m.basic_model_mip, λX, λ, idx_b, idx_a)
            mcbi_post_λxλ_mc(m.basic_model_mip, λλ, λ, idx_a, idx_b)
            mcbi_post_XX_mc(m.basic_model_mip, idx_ab, λX, λλ, lb, ub, idx_a, idx_b)
        end

        # Error condition
        if !(idx_a in m.discrete_x) && !(idx_b in m.discrete_x)
            error("Error case. At least one term should show up in discretization choices.")
        end
    end
end

function mcbi_post_λX(m::JuMP.Model, λX::Dict, dim::Int, idx_a::Int, idx_b::Int)
    λX[(idx_a,idx_b)] = @variable(m, [1:dim], basename=string("L",idx_a,"X",idx_b))
    return λX
end

function mcbi_post_λ(m::JuMP.Model, λ::Dict, dim::Int, idx::Int, relax=false)
    λ[idx] = @variable(m, [1:dim], Bin, basename=string("L",idx))
    @constraint(m, sum(λ[idx]) == 1) # The SOS-1 Constraints, not all MIP solver has SOS feature
    return λ
end

function mcbi_post_λλ(m::JuMP.Model, λλ::Dict, dim_a::Int, dim_b::Int, idx_a::Int, idx_b::Int)
    λλ[(idx_a,idx_b)] = @variable(m, [1:dim_a, 1:dim_b], Bin, basename=string("L",idx_a,"L",idx_b))
    λλ[(idx_b,idx_a)] = λλ[(idx_a,idx_b)]'
    return λλ
end

"""
    Partial partioning for the last mc term
"""
function mcbi_post_λλ(m::JuMP.Model, λλ::Dict, λ::Dict, idx_a::Int, idx_b::Int)
    λλ[(idx_a,idx_b)] = λ[idx_a]
    λλ[(idx_a,idx_b)] = reshape(λλ[(idx_a,idx_b)], (length(λ[idx_a]), 1))
    λλ[(idx_b,idx_a)] = λλ[(idx_a,idx_b)]'
    return λλ
end

function mcbi_post_XX_mc(m, ab, λX, λλ, LB, UB, a, b)
    @assert length(LB[a]) == length(UB[a])
    @assert length(LB[b]) == length(UB[b])
    dim_A = length(LB[a])
    dim_B = length(LB[b])
    @constraint(m, Variable(m, ab) .>= dot(λX[(a,b)],LB[a]) + dot(λX[(b,a)],LB[b]) - reshape(LB[a], (1, dim_A))*λλ[(a,b)]*reshape(LB[b], (dim_B, 1)))
	@constraint(m, Variable(m, ab) .>= dot(λX[(a,b)],UB[a]) + dot(λX[(b,a)],UB[b]) - reshape(UB[a], (1, dim_A))*λλ[(a,b)]*reshape(UB[b], (dim_B, 1)))
	@constraint(m, Variable(m, ab) .<= dot(λX[(a,b)],LB[a]) + dot(λX[(b,a)],UB[b]) - reshape(LB[a], (1, dim_A))*λλ[(a,b)]*reshape(UB[b], (dim_B, 1)))
	@constraint(m, Variable(m, ab) .<= dot(λX[(a,b)],UB[a]) + dot(λX[(b,a)],LB[b]) - reshape(UB[a], (1, dim_A))*λλ[(a,b)]*reshape(LB[b], (dim_B, 1)))
end

function mcbi_post_λxX_mc(m::JuMP.Model, λX::Dict, λ::Dict, ind_λ::Int, ind_X::Int)

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
end

function mcbi_post_λxλ_mc(m::JuMP.Model, λλ::Dict, λ::Dict, ind_A::Int, ind_B::Int)

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
end

function mcbi_post_basic_vars(m::PODNonlinearModel)
    @variable(m.basic_model_mip, x[i=1:m.num_var_orig+m.lifted_x_cont])
    for i in 1:m.num_var_orig
        setcategory(x[i], m.var_type_orig[i])
        setlowerbound(x[i], m.l_var_orig[i])
        setupperbound(x[i], m.u_var_orig[i])
    end
end

function mcbi_post_lifted_aff_constraints(m::PODNonlinearModel)
    for i in 1:m.num_constr_orig
        @constraint(m.basic_model_mip,
            sum(m.lifted_constr_aff_mip[i][:coefs][j]*Variable(m.basic_model_mip, m.lifted_constr_aff_mip[i][:vars][j].args[2]) for j in 1:m.lifted_constr_aff_mip[i][:cnt]) >= m.lifted_constr_aff_mip[i][:rhs])
    end
end

function mcbi_post_lifted_aff_obj(m::PODNonlinearModel)
    @objective(m.basic_model_mip, m.sense_orig, sum(m.lifted_obj_aff_mip[:coefs][i]*Variable(m.basic_model_mip, m.lifted_obj_aff_mip[:vars][i].args[2]) for i in 1:m.lifted_obj_aff_mip[:cnt]))
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
    status = solve(minvertex)

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
