function add_linking_constraints(m::Optimizer, λ::Dict)
    """
    Add linking constraints between λ[i], λ[j] variables.

    Reference information:
        Jongeun Kim, Jean-Philippe P. Richard, Mohit Tawarmalani, Piecewise Polyhedral Relaxations of Multilinear Optimization, 
        http://www.optimization-online.org/DB_HTML/2022/07/8974.html

    For example, suppose we have λ[i], λ[j], and λ[k] where i=(1,2,3), j=(1,2,4), and k=(1,2,5).
    λ[i] contains all multipliers for the extreme points in the space of (x1,x2,x3). 
    λ[j] contains all multipliers for the extreme points in the space of (x1,x2,x4).
    λ[k] contains all multipliers for the extreme points in the space of (x1,x2,x5).

    Using λ[i], λ[j], or λ[k], we can express multilinear function x1*x2.
    We define a linking variable μ(1,2) that represents the value of x1*x2.
    Linking constraints are
        μ(1,2) == convex combination expression for x1*x2 using λ[i],    
        μ(1,2) == convex combination expression for x1*x2 using λ[j], and   
        μ(1,2) == convex combination expression for x1*x2 using λ[k].    

    Then, these constraints link between λ[i], λ[j], and λ[k] variables.
    """

    # (Step 0) Define informational paramters used in the report in the end of the function.
    t_begin = time()
    num_added_constraints = 0

    # (Step 1) Define link_info::Dict 
    #   - key: linking variable indices
    #   - value: the subset of key(λ) where each element is a superset of the linking variable indices 
    maxdeg = maximum([length(k) for k in keys(λ)]) # Compute maximum degree of multilinear terms
    if maxdeg <= 2
        # @info "no linking constraint is added because the maximum degree of multilinear term is $(maxdeg)."
        return
    end
    all_vars_idx = sort(collect(union(collect(keys(λ))...))) # Collect all variable indices used in multilinear terms
    link_info = Dict(
        link_var_idx => filter(r -> issubset(link_var_idx, r), keys(λ)) for
        deg in 2:maxdeg-1 for
        link_var_idx in Combinatorics.combinations(all_vars_idx, deg)
    )
    filter!(r -> length(r.second) >= 2, link_info)
    if isempty(link_info)
        # @info "no linking constraint is added because there is no multilinear terms that shaer more than two variables."
        return
    end

    # (Step 2) Add linking (decision) variables to the mip model
    link_vars = @variable(m.model_mip, [collect(keys(link_info))])

    # (Step 3) Add linking constraints to the mip model
    for (link_var_idx, list_multi_idx) in link_info,
        (i, multi_idx) in enumerate(list_multi_idx)
        # Define `linkidx2multiidx` to know where `link_var_idx` and `list_multi_idx` have the same value.
        # `linkidx2multiidx` satisfies that link_var_idx[i] = list_multi_idx[linkidx2multiidx[i]] for all i. 
        # By the definition of `link_info`, it always holds that issubset(link_var_idx, list_multi_idx).
        linkidx2multiidx = [findfirst(multi_idx .== i) for i in link_var_idx]

        # We define the RHS of linking constraints (convex combination expression for prod(x[link_var_idx]) using λ[multi_idx])
        expr = 0
        for idx in Iterators.product([1:d for d in λ[multi_idx][:dim]]...)
            var = λ[multi_idx][:vars][λ[multi_idx][:indices][idx...]]
            func_value = prod(
                m.discretization[i][idx[j]] for
                (i, j) in zip(link_var_idx, linkidx2multiidx)
            )
            expr += func_value * var
        end

        @constraint(m.model_mip, link_vars[link_var_idx] == expr)
        num_added_constraints += 1
    end

    # @info "$num_added_constraints linking constraints are added in $(time() - t_begin) seconds."
end
