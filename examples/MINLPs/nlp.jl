function nlp1(; solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, 1 <= x[1:2] <= 4)
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

function mathopt(; solver = nothing)
    # Modified based on `mathopt` instances from MINLPLib
    # Global solution: -1004.73094034, arg min: [ -10.0, 15.2543377, -3.2034701]

    m = JuMP.Model(solver)

    LB = [-10, -15, -10]
    UB = [20, 20, 20]
    @variable(m, LB[i] <= x[i = 1:3] <= UB[i])
    @NLobjective(m, Min, 10 * (x[1] * x[2] - x[2] * x[3]) + x[1] * x[3])
    @NLconstraint(m, x[1] * x[2] + x[2] * x[3]^2 >= 4)
    @constraint(m, 3 * x[1] + 4 * x[2] + 5 * x[3] <= 15)

    return m
end

function nlp3(; solver = nothing)
    m = JuMP.Model(solver)

    LB = [100, 1000, 1000, 10, 10, 10, 10, 10]
    UB = [10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000]

    @variable(m, LB[i] <= x[i = 1:8] <= UB[i])

    @constraint(m, 0.0025 * (x[4] + x[6]) <= 1)
    @constraint(m, 0.0025 * (x[5] - x[4] + x[7]) <= 1)
    @constraint(m, 0.01(x[8] - x[5]) <= 1)
    @NLconstraint(m, 100 * x[1] - x[1] * x[6] + 833.33252 * x[4] <= 83333.333)
    @NLconstraint(m, x[2] * x[4] - x[2] * x[7] - 1250 * x[4] + 1250 * x[5] <= 0)
    @NLconstraint(m, x[3] * x[5] - x[3] * x[8] - 2500 * x[5] + 1250000 <= 0)

    @objective(m, Min, x[1] + x[2] + x[3])

    return m
end
