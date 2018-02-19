function sincos_partition_injector(m::PODNonlinearModel, var::Int, partvec::Vector, value::Any, ratio::Any, processed::Set{Int})

    var in processed && return

    # Only consider variables that shows up in :sin/:cos terms, and collect all operators
    operators = Set()
    for nlk in keys(m.nonconvex_terms)
        if m.nonconvex_terms[nlk][:nonlinear_type] in [:sin, :cos]
            if var in m.nonconvex_terms[nlk][:var_idxs]
                push!(operators, m.nonconvex_terms[nlk][:nonlinear_type])
            end
        end
    end
    isempty(operators) && return

    m.var_type[var] != :Cont && error("Method current doesn't support discrete variables in trigonometry functions.")

    λCnt = length(partvec)

    if λCnt == 2  # Initialization
        # Pre-segment the domain based on derivative information (limited to sin/cos)
        seg_points = []
        for o in operators
            (partvec[end] - partvec[1] > 100*pi) && error("Bounds on VAR$(var) is too large $((partvec[end] - partvec[1])/pi)PI. Error this run for performance consideration.")
            # seg_points = [seg_points, tri_zero_der(partvec[1], partvec[end], o);]
            seg_points = [seg_points, tri_convexity_flip_points(partvec[1], partvec[end], o);]
        end
        seg_points = sort(unique(seg_points), rev=true) # pre-segment points
        for i in seg_points
            i in partvec || insert!(partvec, 2, i)
        end
        push!(processed, var)
    else
        point = correct_point(m, partvec, value)
        for j in 1:λCnt
            if point >= partvec[j] && point <= partvec[j+1]  # Locating the right location
                radius = calculate_radius(partvec, j, ratio)
                insert_partition(m, var, j, point, radius, partvec)
                push!(processed, var)
                break
            end
        end
    end

    return
end

"""
    Relative angle location in a 2π range
"""
function tri_loc(x)
    return mod(x, 2*pi)/pi
end

"""
    Derivaive value of trigonometric functions
"""
function tri_der(x, opt)
    opt == :sin && return cos(x)
    opt == :cos && return -sin(x)
    opt == :tan && return sec(x)
    erro("unsupported trigonometric function for derivative")
end

"""
    Trigonometric values
"""
function tri_val(x, opt)
    return eval(opt(x))
end

"""
    Return a vector of values within [a->b] such that the derivative of trigonometric function is 0
"""
function tri_zero_der(a, b, opt)
    opt == :tan && return []

    vals = []

    if opt == :sin
        a_pos = tri_loc(a)
        if a_pos < 0.5
            ex = (0.5-a_pos)*pi
        elseif a_pos < 1.5
            ex = (1.5-a_pos)*pi
        elseif a_pos in [0.5, 1.5]
            ex = 0.0
        else
            ex = (0.5+(2.0-a_pos))*pi
        end
        a + ex > b && return vals
        a + ex < b && push!(vals, a+ex) # Collect the initial value
        nx = a+ex
        while nx <= b
            nx += pi
            nx < b && push!(vals, nx)
        end
        return vals
    elseif opt == :cos
        a_pos = tri_loc(a)
        if a_pos < 1.0 && a_pos > 0.0
            ex = (1.0-a_pos)*pi
        elseif a_pos < 2.0
            ex = (2.0-a_pos)*pi
        elseif a_pos in [0.0, 1.0, 2.0]
            ex = 0.0
        else
            error("EXCEPTION unexpected pos condition in tri_zero_der")
        end
        a + ex > b && return vals
        a + ex < b && push!(vals, a + ex)
        nx = a + ex
        while nx <= b
            nx += pi
            nx < b && push!(vals, nx)
        end
        return vals
    else
        error("Function tri_zero_der currently only supports :sin, :cos, and :tan")
    end

    return vals
end

function tri_convexity_flip_points(a, b, opt)

    opt in [:sin, :cos] || return []

    vals = []

    if opt == :sin
        a_pos = tri_loc(a)
        if a_pos < 1.0
            ex = (1.0-a_pos)*pi
        elseif a_pos in [1.0]
            ex = 0.0
        else
            ex = (2.0-a_pos)*pi
        end
        a + ex > b && return vals
        a + ex < b && push!(vals, a+ex) # Collect the initial value
        nx = a + ex
        while nx <= b
            nx += pi
            nx < b && push!(vals, nx)
        end
        return vals
    elseif opt == :cos
        a_pos = tri_loc(a)
        if a_pos < 0.5
            ex = (0.5-a_pos)*pi
        elseif a_pos < 1.5
            ex = (1.5-a_pos)*pi
        elseif a_pos in [0.5, 1.5]
            ex = 0.0
        else
            ex = (0.5+(2.0-a_pos))*pi
        end
        a + ex > b && return vals
        a + ex < b && push!(vals, a + ex)
        nx = a + ex
        while nx <= b
            nx += pi
            nx < b && push!(vals, nx)
        end
        return vals
    else
        error("Function tri_zero_der currently only supports :sin, :cos, and :tan")
    end

    return vals
end

"""
    Calculate extreme values within range [a, b] and return then by sequence of [a->b]
"""
function tri_extreme_val(a, b, opt)

    if b - a >= 2*pi
        tri_der(a, opt) > 0.0 && return 1, -1
        tri_der(a, opt) < 0.0 && return 1, -1
        tri_der(a, opt) == 0.0 && return tri_val(a), -tri_val(a)
    end

    if tri_der(a, opt) * tri_der(b, opt) > 0.0
        if (b - a) >= pi
            tri_der(a, opt) > 0 && return 1, -1
        else
            tri_der(a, opt) < 0 && return -1, 1
        end
    elseif tri_der(a, opt) * tri_der(b, opt) < 0.0
        tri_der(a, opt) > 0.0 && return 1, min(tri_val(a), tri_val(b))
        tri_der(a, opt) < 0.0 && return -1, max(tri_val(a), tri_val(b))
    else
        if tri_der(a, opt) == 0.0 && tri_der(b, opt) == 0.0

        elseif tri_der(a, opt) == 0.0 && tri_der(b, opt) != 0.0

        elseif tri_der(a, opt) != 0.0 && tri_der(b, opt) == 0.0

        end
    end

end
