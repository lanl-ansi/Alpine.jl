#=====================----- Developer's note's -----=====================#
#  Code used to use convexhull representation to convexify the non-convex model.
#========================================================================#
function amp_post_convhull(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_discretization) ? discretization = options[:use_discretization] : discretization = m.discretization

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables
    Y = Dict()  # Additional variables

    # Construct λ variable space
    for bi in keys(m.nonlinear_terms)
        nl_type = m.nonlinear_terms[bi][:nonlinear_type]
        if ((nl_type == :multilinear) || (nl_type == :bilinear)) && (m.nonlinear_terms[bi][:convexified] == false)
            m.nonlinear_terms[bi][:convexified] = true  # Bookeeping the examined terms
            ml_indices, dim, extreme_point_cnt = amp_convhull_prepare(discretization, bi)   # convert key to easy read mode
            λ = amp_convhull_λ(m, bi, ml_indices, λ, extreme_point_cnt, dim)
            λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, ones(Int,length(dim)))
            α, Y = amp_convhull_α(m, ml_indices, α, Y, dim, discretization)
            amp_post_convhull_constrs(m, λ, α, Y, ml_indices, dim, extreme_point_cnt, discretization)
        elseif (nl_type == :monomial) && (m.nonlinear_terms[bi][:convexified] == false)
            m.nonlinear_terms[bi][:convexified] = true
            monomial_index, dim, extreme_point_cnt = amp_convhull_prepare(discretization, bi, monomial=true)
            λ = amp_convhull_λ(m, bi, monomial_index, λ, extreme_point_cnt, dim)
            λ = populate_convhull_extreme_values(m, discretization, monomial_index, λ)
            α, Y = amp_convhull_α(m, [monomial_index], α, Y, dim, discretization)
            amp_post_convhull_constrs(m, λ, α, Y, monomial_index, dim, discretization)
        end
    end

    return
end

"""
    TODO: docstring
    This function is very important.
"""
function amp_convhull_prepare(discretization::Dict, nonlinear_key::Any; monomial=false)

    id = Set()                      # Coverting the nonlinear indices into a different space
    for var in nonlinear_key        # This output regulates the sequence of how composing variable should be arranged
        push!(id, var.args[2])
    end

    if length(id) < length(nonlinear_key) # Got repeating terms, now the sequence matters
        id = []
        for var in nonlinear_key
            push!(id, var.args[2])
        end
    end

    dim = []
    for i in id
        push!(dim, length(discretization[i]))
    end

    monomial && return id[1], tuple(dim[1]), dim[1]   # One less dimension is required

    return id, tuple([i for i in dim]...), prod(dim)
end

function amp_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, ml_indices::Any, λ::Dict, extreme_point_cnt::Int, dim::Tuple; kwargs...)

    lifted_var_idx = m.nonlinear_terms[nonlinear_key][:lifted_var_ref].args[2]
    @assert !(lifted_var_idx in keys(λ))
    λ[ml_indices] = Dict(:dim=>dim,
                         :lifted_var_idx=>lifted_var_idx,
                         :indices=>reshape([1:extreme_point_cnt;], dim),
                         :vars=>@variable(m.model_mip, [1:extreme_point_cnt], lowerbound=0, upperbound=1, basename="L$(lifted_var_idx)"),
                         :vals=>ones(dim))

    return λ
end

"""
    TODO: docstring
    method for monomial
"""
function populate_convhull_extreme_values(m::PODNonlinearModel, discretization::Dict, monomial_index::Int, λ::Dict)
    λ[monomial_index][:vals] = [discretization[monomial_index][i]^2 for i in 1:length(discretization[monomial_index])]
    return λ
end

"""
    TODO: docstring
    method for general multilinear term
"""
function populate_convhull_extreme_values(m::PODNonlinearModel, discretization::Dict, ml_indices::Any, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

    if level > length(dim)
        @assert length(ml_indices) == length(dim)
        @assert length(ml_indices) == length(locator)
        val = 1.0
        k = 0
        for i in ml_indices
            k += 1
            val *= discretization[i][locator[k]]                      # Calculate extreme point z-value
        end
        λ[ml_indices][:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i
            λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, locator, level+1)
        end
    end

    return λ
end

"""
    Less memeory & time efficient, but a easier implementation
"""
function _populate_convhull_extreme_values(discretization::Dict, ml_indices::Any, λ::Dict, extreme_point_cnt::Int)

    var_indices = collect(ml_indices)
    for i in 1:extreme_point_cnt
        sub = ind2sub(λ[ml_indices][:indices], i)
        k = 0
        for var in ml_indices
            k += 1
            @assert var == var_indices[k]
            λ[ml_indices][:vals][CartesianIndex(sub)] *= discretization[var][sub[k]]
        end
    end

    return λ
end

function amp_convhull_α(m::PODNonlinearModel, ml_indices::Any, α::Dict, Y::Dict, dim::Tuple, discretization::Dict; kwargs...)

    for i in ml_indices
        if !(i in keys(α))
            if m.convhull_formulation_sos2aux
                lambda_cnt = length(discretization[i])
                α[i] = @variable(m.model_mip, [1:intersect_cnt], lowerbound=0.0, upperbound=1.0, basename="A$(i)")
                addSOS2(m.model_mip, α[i])
                @constraint(m.model_mip, Variable(m.model_mip, i) >= sum(α[i][j]*discretization[i][j] for j in 1:lambda_cnt-1)) # Add x = f(α) for regulating the domains
                @constraint(m.model_mip, Variable(m.model_mip, i) <= sum(α[i][j-1]*discretization[i][j] for j in 2:lambda_cnt))
            else
                lambda_cnt = length(discretization[i])
                partition_cnt = length(discretization[i]) - 1
                if m.embedding_sos1
                    YCnt = Int(ceil(log(2,partition_cnt)))
                    α[i] = @variable(m.model_mip, [1:partition_cnt], lowerbound=0.0, upperbound=1.0, basename="A$(i)")
                    Y[i] = @variable(m.model_mip, [1:YCnt], Bin, basename=string("YL",i))
                    αYmap = ebd_sos1(partition_cnt, m.embedding_encode)
                    for j in 1:YCnt  # Setup the constraint to link α and Y
                        @constraint(m.model_mip, sum(α[i][k] for k=1:partition_cnt if k in αYmap[j]) <= Y[i][Int(j)])
                        @constraint(m.model_mip, sum(α[i][k] for k=1:partition_cnt if k in αYmap[YCnt+j]) <= 1 - Y[i][Int(j)])
                    end
                    @constraint(m.model_mip, sum(α[i]) == 1)
                    @constraint(m.model_mip, Variable(m.model_mip, i) >= sum(α[i][j]*discretization[i][j] for j in 1:lambda_cnt-1)) # Add x = f(α) for regulating the domains
                    @constraint(m.model_mip, Variable(m.model_mip, i) <= sum(α[i][j-1]*discretization[i][j] for j in 2:lambda_cnt))
                elseif m.embedding_sos2
                    αCnt = Int(ceil(log(2,partition_cnt)))
                    α[i] = @variable(m.model_mip, [1:αCnt], Bin, basename=string("YL",i))
                else
                    α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(i)")
                    @constraint(m.model_mip, sum(α[i]) == 1)
                    @constraint(m.model_mip, Variable(m.model_mip, i) >= sum(α[i][j]*discretization[i][j] for j in 1:lambda_cnt-1)) # Add x = f(α) for regulating the domains
                    @constraint(m.model_mip, Variable(m.model_mip, i) <= sum(α[i][j-1]*discretization[i][j] for j in 2:lambda_cnt))
                end
            end
        end
    end

    return α, Y
end

function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, Y::Dict, ml_indices::Any, dim::Tuple, extreme_point_cnt::Int, discretization::Dict)

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[ml_indices][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[ml_indices][:lifted_var_idx]) == dot(λ[ml_indices][:vars], reshape(λ[ml_indices][:vals], extreme_point_cnt)))

    cnt = 0
    for i in ml_indices
        cnt += 1
        valid_inequalities(m, discretization, λ, α, Y, ml_indices, dim, i, cnt)        # Add links between λ and α
        if m.convhull_formulation_sos2aux
            @constraint(m.model_mip, Variable(m.model_mip, i) == dot(α[i], discretization[i]))
        else
            partition_cnt = length(α[i])
            lambda_cnt = length(discretization[i])
            @assert lambda_cnt == partition_cnt + 1
            sliced_indices = [collect_indices(λ[ml_indices][:indices], cnt, [k], dim) for k in 1:lambda_cnt] # Add x = f(λ) for convex representation of x value
            @constraint(m.model_mip, Variable(m.model_mip, i) == sum(dot(repmat([discretization[i][k]],length(sliced_indices[k])), λ[ml_indices][:vars][sliced_indices[k]]) for k in 1:lambda_cnt))
        end
    end

    return
end

function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, Y::Dict, monomial_idx::Int, dim::Tuple, discretization::Dict)

    partition_cnt = length(α[monomial_idx])
    lambda_cnt = length(discretization[monomial_idx])
    @assert lambda_cnt == partition_cnt + 1

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[monomial_idx][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[monomial_idx][:lifted_var_idx]) <= dot(λ[monomial_idx][:vars], λ[monomial_idx][:vals]))
    @constraint(m.model_mip, Variable(m.model_mip, λ[monomial_idx][:lifted_var_idx]) >= Variable(m.model_mip, monomial_idx)^2)

    # Add SOS-2 Constraints with basic encoding
    for i in 1:lambda_cnt
        if i == 1
            @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i])
        elseif i == lambda_cnt
            @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1])
        else
            @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1] + α[monomial_idx][i])
        end
    end

    # Equalivent SOS-2 : A different encoding (solwer performance)
    # for i in 1:partition_cnt
    #     @constraint(m.model_mip, α[monomial_idx][i] <= λ[monomial_idx][:vars][i] + λ[monomial_idx][:vars][i+1])
    # end
    # @constraint(m.model_mip, α[monomial_idx][1] >= λ[monomial_idx][:vars][1])
    # @constraint(m.model_mip, α[monomial_idx][end] >= λ[monomial_idx][:vars][end])

    # Add x = f(λ) for convex representation
    @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) == dot(λ[monomial_idx][:vars], discretization[monomial_idx]))

    # Add x = f(α) for regulating the domains
    @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) >= sum(α[monomial_idx][j]*discretization[monomial_idx][j] for j in 1:lambda_cnt-1))
    @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) <= sum(α[monomial_idx][j-1]*discretization[monomial_idx][j] for j in 2:lambda_cnt))

    (m.convhull_formulation_sos2aux) && warn("Not considering alternative SOS2-Formulation for monomial terms")

    return
end

# Valid inequalities proposed when Jeff L. was here
function valid_inequalities(m::PODNonlinearModel, discretization::Dict, λ::Dict, α::Dict, Y::Dict, ml_indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    partition_cnt = length(α[var_ind])
    lambda_cnt = length(discretization[var_ind])

    if m.convhull_formulation_sos2aux
        for j in 1:lambda_cnt
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j], dim)
            @constraint(m.model_mip, α[var_ind][j] == sum(λ[ml_indices][:vars][sliced_indices]))
        end
        return
    end

    if m.convhull_formulation_facet
        # Constraint cluster of α >= f(λ)
        for j in 1:(partition_cnt-1) # Construct cuts by sweeping in both directions
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        end

        # Constriant cluster of α <= f(λ)
        for j in 1:(partition_cnt-1)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:(j+1);], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) <= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        return
    end

    if m.convhull_formulation_sos2
        # Encoding of λ -> α goes here
        if m.embedding_sos2
            ebd_map = ebd_sos2(lambda_cnt, m.embedding_encode)
            YCnt = length(keys(ebd_map))/2
            for i in 1:YCnt
                p_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, ebd_map[var_ind], dim)
                n_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, ebd_map[var_ind+YCnt], dim)
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][p_sliced_indices]) <= α[var_ind][i])
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][n_sliced_indices]) <= 1-α[var_ind][i])
            end
        else
            for j in 1:lambda_cnt
                sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j], dim)
                if (j == 1)
                    @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][j])
                elseif (j == lambda_cnt)
                    @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][partition_cnt])
                else
                    @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= sum(α[var_ind][(j-1):j]))
                end
            end
        end
        return
    end

    if m.convhull_formulation_minib

        # Constraint cluster of α >= f(λ)
    	for j in 1:min(partition_cnt, m.convexhull_sweep_limit) # Construct cuts by sweeping in both directions
        	sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
        	@constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        	sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [(lambda_cnt-j+1):(lambda_cnt);], dim)
        	@constraint(m.model_mip, sum(α[var_ind][(dim[cnt]-j):(dim[cnt]-1)]) >= sum(λ[ml_indices][:vars][sliced_indices]))
    	end

    	# Constriant cluster of α <= f(λ)
    	for j in 1:partition_cnt
       		for i in 1:max(1, min(partition_cnt-j+1, m.convexhull_sweep_limit)) # At least one
            	sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j:(j+i);], dim)
            	@constraint(m.model_mip, sum(α[var_ind][j:(j+i-1)]) <= sum(λ[ml_indices][:vars][sliced_indices]))
        	end
    	end

        return
    end

    error("Must indicate a choice of convex hull formulation. ?(minib, sos2, sos2aux, facet)")
    return
end

function collect_indices(l::Array, locator::Tuple, dim::Tuple)

    k = 0
    indices = Vector{Int}(2^length(dim))
    for i in 1:prod(dim)
        ind = ind2sub(l, i)
        diff = [((ind[i] - locator[i] == 0) || (ind[i] - locator[i] == 1)) for i in 1:length(dim)]
        if prod(diff)
            k +=1
            indices[k] = i
        end
    end

    return indices
end

function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Array, dim::Tuple)

	k = 0
	indices = Vector{Int}(Int(prod(dim)/dim[fixed_dim]*length(fixed_partition)))
	for i in 1:prod(dim)
		ind = ind2sub(l, i)
		if ind[fixed_dim] in fixed_partition
			k += 1
			indices[k] = i
		end
	end

	return indices
end

function resolve_lifted_var_value(m::PODNonlinearModel, sol_vec::Array)

    @assert length(sol_vec) == m.num_var_orig
    sol_vec = [sol_vec; fill(NaN, m.num_var_linear_lifted_mip+m.num_var_nonlinear_lifted_mip)]

    for i in 1:length(m.nonlinear_terms)
        for bi in keys(m.nonlinear_terms)
            if m.nonlinear_terms[bi][:id] == i
                lvar_idx = m.num_var_orig + i
                if haskey(m.nonlinear_terms[bi], :evaluator)
                    sol_vec[lvar_idx] = m.nonlinear_terms[bi][:evaluator](m.nonlinear_terms[bi], sol_vec)
                else
                    sol_vec[lvar_idx] = 0.5*m.discretization[lvar_idx][1] + 0.5*m.discretization[lvar_idx][end]
                end
            end
        end
    end

    return sol_vec
end
