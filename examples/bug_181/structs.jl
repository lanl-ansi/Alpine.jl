abstract type Data end

struct TranspProblemData <: Data
    caps::Array{Int64, 1}
    demands::Array{Int64, 1}
    asscosts::Array{Int64, 2}
end
