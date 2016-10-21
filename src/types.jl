import Base: (*), (+), (-), (/), (^), (==)

abstract Xial

type Monomial <: Xial
    c::Float64
    terms::Dict{Int,Float64} # variable index to power (coefficient)
end
Monomial(c::Number) = Monomial(c, Dict{Int,Float64}())
Base.convert(::Type{Monomial}, c::Number) = Monomial(c)
(==)(m1::Monomial, m2::Monomial) =
(m1.c == m2.c) && (m1.terms == m2.terms)
function Base.print(io::IO, mon::Monomial)
    print(io, mon.c)
    for (i,v) in mon.terms
        print(io, "*x_{$i}^{$v}")
    end
end


type Polynomial <: Xial
    mons::Vector{Monomial}
end

Polynomial(c::Number) = Polynomial(Monomial(c))
Polynomial(m::Monomial) = Polynomial(Monomial[m])
Polynomial(args...) = Polynomial(Monomial[args...])
function (==)(p1::Polynomial, p2::Polynomial)
    # Doesn't handle sorting differences
    length(p1.mons) != length(p2.mons) && return false
    all(pair->pair[1]==pair[2], zip(p1.mons,p2.mons))
end
function Base.print(io::IO, pos::Polynomial)
    if length(pos.mons) == 0
        print(io, "0")
    elseif length(pos.mons) == 1
        print(io, pos.mons[1])
    else
        print(io, "[")
        for mon in pos.mons[1:end-1]
            print(io, mon, " + ")
        end
        print(io, pos.mons[end], "]")
    end
end

# other types go here
