"""
    Find an integer point from solution point. Wrapps discover_else_integer() when rounded solution
    step on a integer that has already been discovered.
"""
function discover_int_point(m::PODNonlinearModel, var::Int, partvec::Vector, point::Float64)

    if point < partvec[1] + m.tol || point > partvec[end] - m.tol
        warn("Integer Soluiton SOL=$(point) out of bounds [$(d[i][1]),$(d[i][end])].")
        return discover_else_integer(m, var, partvec)
    end

    round_int = Int(round(point))

    # Only return the value that is still convered within relaxed region
    if is_uncovered_integer(round_int, partvec)
        return discover_else_integer(m, var, partvec)
    else
        return round_int
    end
end

"""
    Find an integer point from the largest uncovered discretized region
    If variable domain is fully discovered, return nothing
"""
function discover_else_integer(m::PODNonlinearModel, var::Int, partvec::Vector)

    l_cnt = length(partvec)
    p_cnt = length(partvec) - 1

    point = nothing

    incumblength = 0.0
    incumbpart = -1
    for i in 1:p_cnt
        if partvec[i+1] - partvec[i] > incumblength
            incumblength = partvec[i+1] - partvec[i]
            incumbpart = i
        end
    end

    # Region has been fully discovered
    incumblength == 1.0 && return nothing
    interval_cnt = (partvec[incumbpart+1]) - (partvec[incumbpart])

    if mod(interval_cnt, 2) == 0
        return partvec[incumbpart] + 0.5 + interval_cnt/2 - 1
    else
        return partvec[incumbpart] + interval_cnt/2
    end

    m.loglevel > 99 && prinln("[EXIT] NEW POINT = $(point)")
    return point
end

function is_uncovered_integer(val::Int, partvec::Vector)

    if val - 0.5 in partvec && val + 0.5 in partvec
        return true
    else
        return false
    end
end
