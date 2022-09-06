"""
    relaxation_bilinear(
    m::JuMP.Model, 
    z::JuMP.VariableRef, 
    x::JuMP.VariableRef, 
    y::JuMP.VariableRef, 
    lb_x::Number,
    ub_x::Number,
    lb_y::Number,
    ub_y::Number)

General relaxation of a binlinear term (McCormick relaxation). 
```
z >= JuMP.lower_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.lower_bound(x)*JuMP.lower_bound(y)
z >= JuMP.upper_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.upper_bound(x)*JuMP.upper_bound(y)
z <= JuMP.lower_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.lower_bound(x)*JuMP.upper_bound(y)
z <= JuMP.upper_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.upper_bound(x)*JuMP.lower_bound(y)
```
"""
function relaxation_bilinear(
    m::JuMP.Model,
    z::JuMP.VariableRef,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    lb_x::Number,
    ub_x::Number,
    lb_y::Number,
    ub_y::Number
)

    JuMP.@constraint(m, z >= lb_x * y + lb_y * x - lb_x * lb_y)
    JuMP.@constraint(m, z >= ub_x * y + ub_y * x - ub_x * ub_y)
    JuMP.@constraint(m, z <= lb_x * y + ub_y * x - lb_x * ub_y)
    JuMP.@constraint(m, z <= ub_x * y + lb_y * x - ub_x * lb_y)

    return
end

"""
    relaxation_multilinear_binary(
    m::JuMP.Model, 
    z::JuMP.VariableRef, 
    x::Vector{VariableRef})

Applies Fortet linearization (see https://doi.org/10.1007/s10288-006-0015-3) for z = prod(x), 
where x is a vector of binary variables.
"""
function relaxation_multilinear_binary(m::JuMP.Model, z::JuMP.VariableRef, x::Vector{VariableRef})
    for i in x
        JuMP.@constraint(m, z <= i)
    end
    JuMP.@constraint(m, z >= sum(x) - (length(x) - 1))

    return
end

#=
"""
    relaxation_quadratic_univariate(
        m::JuMP.Model,
        z::JuMP.VariableRef,
        x::JuMP.VariableRef,
        lb_x::Number,
        ub_x::Number,
    )

Applies convex hull relaxation for z = x^2, where x is a variable.
"""
function relaxation_quadratic_univariate(
    m::JuMP.Model,
    z::JuMP.VariableRef,
    x::JuMP.VariableRef,
    lb_x::Number,
    ub_x::Number,
)
    JuMP.@constraint(m, z >= x^2)
    JuMP.@constraint(m, z <= (lb_x + ub_x) * x - (lb_x * ub_x))

    return
end
=#