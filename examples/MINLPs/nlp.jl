function nlp1(; solver = nothing)
    m = JuMP.Model(solver)

    sol = [2.55577, 3.13017]
    @variable(m, 1 <= x[i = 1:2] <= 4, start = sol[i])
    @NLconstraint(m, x[1] * x[2] >= 8)
    @NLobjective(m, Min, 6 * x[1]^2 + 4 * x[2]^2 - 2.5 * x[1] * x[2])

    return m
end

function nlp2(; solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, -500 <= x[1:2] <= 500)
    @NLobjective(m, Min, sum((x[i]^2 - i)^2 for i in 1:2))

    return m
end

function max_cover_var_picker(m::Alpine.Optimizer)
    nodes = Set()
    for pair in keys(m.nonconvex_terms)
        for i in pair
            @assert isa(i.args[2], Int)
            push!(nodes, i.args[2])
        end
    end
    nodes = collect(nodes)
    m.num_var_disc_mip = length(nodes)
    m.disc_vars = nodes
    return
end

function nlp3(; solver = nothing)
    m = JuMP.Model(solver)

    LB = [100, 1000, 1000, 10, 10, 10, 10, 10]
    UB = [10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000]

    sol = [579.30669, 1359.97066, 5109.97054, 182.0177, 295.60118, 217.9823, 286.41653, 395.60118]
    @variable(m, LB[i] <= x[i = 1:8] <= UB[i], start = sol[i])

    @constraint(m, 0.0025 * (x[4] + x[6]) <= 1)
    @constraint(m, 0.0025 * (x[5] - x[4] + x[7]) <= 1)
    @constraint(m, 0.01(x[8] - x[5]) <= 1)
    @NLconstraint(m, 100 * x[1] - x[1] * x[6] + 833.33252 * x[4] <= 83333.333)
    @NLconstraint(m, x[2] * x[4] - x[2] * x[7] - 1250 * x[4] + 1250 * x[5] <= 0)
    @NLconstraint(m, x[3] * x[5] - x[3] * x[8] - 2500 * x[5] + 1250000 <= 0)

    @objective(m, Min, x[1] + x[2] + x[3])

    return m
end
