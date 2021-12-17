using NLsolve

function cubic_function_x_multipliers(maxx)
    return maxx * [0, 1/3, 2/3, 1]
end

function cubic_function_y_multipliers(nomcost)
    return nomcost * [1/2, 3/4, 7/8, 2]
end

function cubic_function_multipliers(maxx, nomcost)
    cubic_function_x_multipliers(maxx), cubic_function_y_multipliers(nomcost)
end

function cubic_function(x, q)
    return x[1] * q^3 + x[2] * q^2 + x[3] * q + x[4]
end

function cubic_coefs(x::Array{Float64, 1}, y::Array{Float64, 1})
    function f!(F, w)
        for i in 1 : 4
            F[i] = cubic_function(w, x[i]) - y[i]
        end
    end
    zeros = nlsolve(f!, [1.0, 1.0, 1.0, 1.0]).zero
end
