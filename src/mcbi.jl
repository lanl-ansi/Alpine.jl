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

function mcbi_pick_discretize_vars(m::PODNonlinearModel)
    # Figure out which are the variables that needs to be partitioned
    if m.discrete_vars_choice == 0
        mcbi_max_cover(m)
    elseif m.discrete_vars_choice == 1
        mcbi_min_vertex_cover(m)
    else
        error("Unsupported method for picking variables for discretization")
    end
end

function mcbi_post_mc_all(m::PODNonlinearModel; kwargs...)

    λ = Dict()
    λλ = Dict()
    λX = Dict()

    for bi in keys(m.dict_nonlinear_info)
        @assert length(bi) == 2
        idx_a = bi[1].args[2]
        idx_b = bi[2].args[2]
        idx_ab = m.dict_nonlinear_info[bi].args[2]
        if (idx_a in m.discrete_x) && !(idx_b in m.discrete_x)
            λ[idx_a] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_a])-1)], Bin, basename=string("l",idx_a))
            λX[(idx_a,idx_b)] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_a])-1)], basename=string("l",idx_a,"x",idx_b))
            λX[(idx_b,idx_a)] = Variable(m.basic_model_mip, idx_b)
            λλ[(idx_a,idx_b)] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_a])-1), 1], Bin, basename=string("l",idx_a,"l",idx_b))
            mcbi_λxX_mc(m.basic_model_mip, λX, λ, X, idx_a, idx_b)
            mcbi_XX_mc(m.basic_model_mip, idx_ab, λX, λλ, ?, ?, idx_a, idx_b)
        end
        if !(idx_a in m.discrete_x) && (idx_b in m.discrete_x)
            λ[idx_b] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_b])-1)], Bin, basename=string("l",idx_b))
            λX[(idx_b,idx_a)] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_b])-1)], basename=string("l",idx_b,"x",idx_a))
            λX[(idx_a,idx_b)] = Variable(m.basic_model_mip, idx_a)
            λλ[(idx_a,idx_b)] = @variable(m.basic_model_mip, [1, 1:(length(m.discretization[idx_b])-1)], Bin, basename=string("l",idx_a,"l",idx_b))
            mcbi_λxX_mc(m, λX, λ, X, idx_b, idx_a)
            mcbi_XX_mc(m.basic_model_mip,idx_ab, λX, λλ, ?, ?, idx_a, ix_b)
        end
        if (idx_a in m.discrete_x) && (idx_b in m.discrete_x)
            λ[idx_a] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_a])-1)], Bin, basename=string("l",idx_a))
            λ[idx_b] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_b])-1)], Bin, basename=string("l",idx_a))
            λX[(idx_a,idx_b)] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_a])-1)], basename=string("l",idx_a,"x",idx_b))
            λX[(idx_b,idx_a)] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_b])-1)], basename=string("l",idx_b,"x",idx_a))
            λλ[(idx_a,idx_b)] = @variable(m.basic_model_mip, [1:(length(m.discretization[idx_a])-1), 1:(length(m.discretization[idx_b])-1)], Bin, basename=string("l",idx_a,"l",idx_b))
            mcbi_λxX_mc(m, λX, λ, X, idx_a, idx_b)
            mcbi_λxX_mc(m, λX, λ, X, idx_b, idx_a)
            mcbi_XX_mc(m.basic_model_mip, idx_ab, λX, λλ, ?, ?, idx_a, ix_b)
        end
        if !(idx_a in m.discrete_x) && !(idx_b in m.discrete_x)
            error("Error case. At least one term should show up in discretization choices.")
        end
    end
end

function mcbi_XX_mc(m, ab, λX, λλ, LB, UB, a, b)
	@constraint(m, Variable(m, ab) >= λX[(a,b)] * LB[a] + λX[(a,b)] * LB[b] - LB[a]' * λλ[(a,b)] * LB[b])
	@constraint(m, Variable(m, ab) >= λX[(a,b)] * UB[a] + λX[(a,b)] * UB[b] - UB[a]' * λλ[(a,b)] * UB[b])
	@constraint(m, Variable(m, ab) <= λX[(a,b)] * LB[a] + λX[(a,b)] * UB[b] - LB[a]' * λλ[(a,b)] * UB[b])
	@constraint(m, Variable(m, ab) <= λX[(a,b)] * UB[a] + λX[(a,b)] * LB[b] - UB[a]' * λλ[(a,b)] * LB[b])
	return m
end

function mcbi_λxX_mc(m::JuMP.Model, λX::Dict, λ::Dict, X::Dict, ind_λ::Int, ind_X::Int)

    # X_u and λ here are vectors, and X is one variable,
	# After the product, we have a polynomial to multiply X,
	# forming |λ| new variables stored in λX dictionary

	dim_λ = length(λ[ind_λ]); # This is how many new variables to be generated
	# Now construct McCormick constraint for these new variables
	for i in 1:dim_λ
		lb_X = getlowerbound(Variable(m, ind_X))
		ub_X = getupperbound(Variable(m, ind_X))
		lb_λ = getlowerbound(λ[ind_λ][i])
		ub_λ = getupperbound(λ[ind_λ][i])
        @assert (lb_λ == 0.0) && (ub_λ == 1.0)
		mccormick(m, λX[(ind_λ,ind_X)][i], λ[ind_λ][i], X[ind_X], lb_λ, ub_λ, lb_X, ub_X)
	end

end

function mcbi_λxλ_mc(m::JuMP.Model, λλ::Dict, λ::Dict, ind_A::Int, ind_B::Int)

	#=  General Case :: A 3 x 2 example, some new variables are formed
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
			mccormick_bin(m, λλ[(ind_A,ind_B)[i,j], λ[ind_A][i], λ[ind_B][j])
		end
	end
	return m, λλ
end


"""
    Perform a minimum vertex cover for selecting variables for discretization
"""
function mcbi_min_vertex_cover(m::PODNonlinearModel; kwargs...)

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
    m.discrete_x_cnt = Int(sum(xVal))
    m.discrete_x = [i for i in nodes if xVal[i]>0]

end

"""
    Collect all envolved variables as a max cover for generate discretization
"""
function mcbi_max_cover(m::PODNonlinearModel; kwargs...)
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
