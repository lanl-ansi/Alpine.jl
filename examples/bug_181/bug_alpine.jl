using JuMP, Alpine, BARON

include("structs.jl")
include("costs.jl")
include("optimizers.jl")
include("reader.jl")

function solve(filename)
    pdata = read_problem(filename)
    
    model = JuMP.Model(optimizer_with_attributes(Alpine.Optimizer, 
                    "nlp_solver" => get_ipopt(),
                    "mip_solver" => get_cplex(),
                    "minlp_solver" => get_juniper(),
                    "time_limit" => 3600,
                    "presolve_bt" => true,
                    "presolve_bt_max_iter" => 7)
                    )

    # model = JuMP.Model(optimizer_with_attributes(BARON.Optimizer, 
    #                                    "maxtime" => 3600, 
    #                                    "CplexLibName" => "/Users/harsha/Documents/CPLEX_20_1/cplex/bin/x86-64_osx/libcplex2010.jnilib"))

    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    @variable(model, 0 <= x[i in 1 : ndeps, j in 1 : ncusts] <= min(pdata.caps[i], pdata.demands[j]))
    @variable(model, y[1 : ndeps, 1 : ncusts], Bin)
    @constraint(model, outflow[i in 1 : ndeps], sum(x[i, j] for j in 1 : ncusts) <= pdata.caps[i])
    @constraint(model, inflow[j in 1 : ncusts], sum(x[i, j] for i in 1 : ndeps) >= pdata.demands[j])
    @constraint(model, link[i in 1 : ndeps, j in 1 : ncusts], x[i, j] - min(pdata.demands[j], pdata.caps[i]) * y[i, j] <= 0)
    dim = 4
    costmult = 1
    cf = zeros(ndeps, ncusts, dim)
    for i in 1 : ndeps, j in 1 : ncusts
        maxoffer = min(pdata.demands[j], pdata.caps[i])
        xf, yf = cubic_function_multipliers(maxoffer, costmult * pdata.asscosts[i, j])
        cf[i, j, :] = cubic_coefs(xf, yf)
    end
    
    @NLexpression(model, obj, sum(cf[i, j, 1] * x[i, j]^3 + cf[i, j, 2] * x[i, j]^2 + cf[i, j, 3] * x[i, j] + cf[i, j, 4] * y[i, j] for i in 1 : ndeps, j in 1 : ncusts))
    @NLobjective(model, Min, obj)
    
    JuMP.optimize!(model)
    JuMP.objective_value(model)
end