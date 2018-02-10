function sincos_p1(;solver=nothing)

    m = Model(solver=solver)

    @variable(m, 0.1<=x[1:10]<=10)
    @NLconstraint(m, [i in 1:9], x[i]*x[i+1]*sin(x[i])>=0.32)
    @NLconstraint(m, [i in 1:9], x[i]*x[i+1]*sin(x[i])<=0.42)
    @NLobjective(m, Min, sum(cos(x[i])+sin(x[i]) for i in 1:10))

    # [osx] Ipopt Result : 8.42226986283713
    # [osx] Ipopt Solution :    0.886654
                             #  0.611244
                             #  0.912239
                             #  0.582146
                             #  0.999768
                             #  0.499317
                             #  1.33843
                             #  0.322466
                             #  3.13137
                             #  10.0

    # [osx] Couenne Results : -10.376627601918297
    # [osx] Couenne Solution : 9.42117
                             # 9.42117
                             # 9.42117
                             # 9.42117
                             # 9.42117
                             # 9.42117
                             # 9.42117
                             # 9.42117
                             # 9.41613
                             # 3.92856

    return m
end
