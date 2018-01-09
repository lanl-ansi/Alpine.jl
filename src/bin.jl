function binprod_relax(m, z, x::Vector)
    for i in x
        @constraint(m, z <= i)
    end
    @constraint(m, z >= sum(x) - (length(x)-1))
    return
end
