# Number--
# --Monomial
(+)(x::Number, mon::Monomial) = Monomial(x) + mon
(-)(x::Number, mon::Monomial) = Monomial(x) - mon
(*)(x::Number, mon::Monomial) = Monomial(x*mon.c, mon.terms)
(/)(x::Number, mon::Monomial) = x*mon^-1
# --Polynomial
(+)(x::Number, pol::Polynomial) = Monomial(x) + pol
(-)(x::Number, pol::Polynomial) = Monomial(x) - pol
(*)(x::Number, pol::Polynomial) = Polynomial(map(mon->(x*mon), pol.mons))
(/)(x::Number, pol::Polynomial) = error("Can't divide a scalar by a Polynomial")

# Monomial--
(+)(mon::Monomial) =  1*mon
(-)(mon::Monomial) = -1*mon
# --Number
(+)(mon::Monomial, x::Number) = mon + Monomial(x)
(-)(mon::Monomial, x::Number) = mon - Monomial(x)
(*)(mon::Monomial, x::Number) = Monomial(mon.c*x, mon.terms)
(/)(mon::Monomial, x::Number) = Monomial(mon.c/x, mon.terms)
(^)(m::Monomial, num::Integer) = Monomial(m.c^num,  # for ambiguity
    Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))
(^)(m::Monomial, num::Number) = Monomial(m.c^num,
    Dict{Int,Float64}([i => m.terms[i]*num for i in keys(m.terms)]))
# --Monomial
(+)(m_1::Monomial, m_2::Monomial) = Polynomial([m_1, m_2])
(-)(m_1::Monomial, m_2::Monomial) = Polynomial([m_1,-m_2])
function (*)(m_1::Monomial, m_2::Monomial)
    d = copy(m_1.terms)
    for (i,v) in m_2.terms
        d[i] = get(d,i,0.0) + v
    end
    return Monomial(m_1.c*m_2.c, d)
end
function (/)(m_1::Monomial, m_2::Monomial)
    d = copy(m_1.terms)
    for (i,v) in m_2.terms
        d[i] = get(d,i,0.0) - v
    end
    return Monomial(m_1.c/m_2.c, d)
end
# -- Polynomial
(+)(mon::Monomial,pol::Polynomial) = Polynomial(vcat(mon,pol.mons))
(-)(mon::Monomial,pol::Polynomial) = Polynomial(vcat(mon,map(-,pol.mons)))
(*)(mon::Monomial,pol::Polynomial) = Polynomial([mon*pm for pm in pol.mons])
(/)(mon::Monomial,pol::Polynomial) = error("Can't divide a Monomial by a Polynomial")

# Polynomial--
(+)(pol::Polynomial) =  1*pol
(-)(pol::Polynomial) = -1*pol
# --Number
(+)(pol::Polynomial, x::Number) = Polynomial(vcat(pol.mons,Monomial( x)))
(-)(pol::Polynomial, x::Number) = Polynomial(vcat(pol.mons,Monomial(-x)))
(*)(pol::Polynomial, x::Number) = Polynomial([pm*x for pm in pol.mons])
(/)(pol::Polynomial, x::Number) = Polynomial([pm/x for pm in pol.mons])
# --Monomial
(+)(pol::Polynomial, mon::Monomial) = Polynomial(vcat(pol.mons, mon))
(-)(pol::Polynomial, mon::Monomial) = Polynomial(vcat(pol.mons,-mon))
(*)(pol::Polynomial, mon::Monomial) = Polynomial([pm*mon for pm in pol.mons])
(/)(pol::Polynomial, mon::Monomial) = Polynomial([pm/mon for pm in pol.mons])
# --Polynomial
(+)(p_1::Polynomial, p_2::Polynomial) = Polynomial(vcat(p_1.mons,p_2.mons))
(-)(p_1::Polynomial, p_2::Polynomial) = Polynomial(vcat(p_1.mons,(-p_2).mons))
function (*)(p_1::Polynomial, p_2::Polynomial)
    mon_out = Monomial[]
    for m_1 in p_1.mons, m_2 in p_2.mons
        push!(mon_out, m_1*m_2)
    end
    return Polynomial(mon_out)
end
(/)(p_1::Polynomial, p_2::Polynomial) = error("Can't divide a Polynomial by a Polynomial")
