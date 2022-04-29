function ex1225a(;solver = nothing)
    # GOpt: 131670.377903
    # This problem may be numerically sensitive 

    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    @variable(m, x[x_Idx] >= 0)
    b_Idx = Any[22, 23, 24]
    @variable(m, 0 <= b[b_Idx] <= 1, Bin)
    i_Idx = Any[16, 17, 18, 19, 20, 21]
    @variable(m, 0 <= i[i_Idx] <= 3, Int)

    JuMP.set_upper_bound(x[1], 80.0)
    JuMP.set_upper_bound(x[2], 25.0)
    JuMP.set_upper_bound(x[3], 45.0)
    JuMP.set_upper_bound(x[4], 2950.0)
    JuMP.set_upper_bound(x[5], 2950.0)
    JuMP.set_upper_bound(x[6], 2950.0)
    JuMP.set_upper_bound(x[7], 400.0)
    JuMP.set_upper_bound(x[8], 400.0)
    JuMP.set_upper_bound(x[9], 400.0)
    JuMP.set_upper_bound(x[10], 350.0)
    JuMP.set_upper_bound(x[11], 350.0)
    JuMP.set_upper_bound(x[12], 350.0)
    JuMP.set_upper_bound(x[13], 1.0)
    JuMP.set_upper_bound(x[14], 1.0)
    JuMP.set_upper_bound(x[15], 1.0)


    # ----- Constraints ----- #
    @NLconstraint(m, e1, -((6329.03+1800*x[1])*i[16]*i[19]*b[22]+(2489.31+1800*x[2])*i[17]*i[20]*b[23]+(3270.27+1800*x[3])*i[18]*i[21]*b[24])+objvar == 0.0)
    @NLconstraint(m, e2, (-19.9* (0.000338983050847458*x[4])^3)-0.161* (0.000338983050847458*x[4])^2*x[10]+1.90169491525424e-7*x[4]* (x[10])^2+x[1] == 0.0)
    @NLconstraint(m, e3, (-1.21* (0.000338983050847458*x[5])^3)-0.0644* (0.000338983050847458*x[5])^2*x[11]+1.91186440677966e-7*x[5]* (x[11])^2+x[2] == 0.0)
    @NLconstraint(m, e4, (-6.52* (0.000338983050847458*x[6])^3)-0.102* (0.000338983050847458*x[6])^2*x[12]+7.86440677966102e-8*x[6]* (x[12])^2+x[3] == 0.0)
    @NLconstraint(m, e5, (-0.00023593220338983*x[4]*x[10])-629* (0.000338983050847458*x[4])^2+0.0116* (x[10])^2+x[7] == 0.0)
    @NLconstraint(m, e6, (-0.001*x[5]*x[11])-215* (0.000338983050847458*x[5])^2+0.115* (x[11])^2+x[8] == 0.0)
    @NLconstraint(m, e7, (-0.000179661016949153*x[6]*x[12])-361* (0.000338983050847458*x[6])^2+0.00946* (x[12])^2+x[9] == 0.0)
    @constraint(m, e8, x[13]+x[14]+x[15] == 1.0)
    @NLconstraint(m, e9, 0.00285714285714286*x[10]*i[16]-x[13] == 0.0)
    @NLconstraint(m, e10, 0.00285714285714286*x[11]*i[17]-x[14] == 0.0)
    @NLconstraint(m, e11, 0.00285714285714286*x[12]*i[18]-x[15] == 0.0)
    @NLconstraint(m, e12, 0.0025*x[7]*i[19]-b[22] == 0.0)
    @NLconstraint(m, e13, 0.0025*x[8]*i[20]-b[23] == 0.0)
    @NLconstraint(m, e14, 0.0025*x[9]*i[21]-b[24] == 0.0)
    @constraint(m, e15, 0.000338983050847458*x[4]-b[22] <= 0.0)
    @constraint(m, e16, 0.000338983050847458*x[5]-b[23] <= 0.0)
    @constraint(m, e17, 0.000338983050847458*x[6]-b[24] <= 0.0)
    @constraint(m, e18, 0.0125*x[1]-b[22] <= 0.0)
    @constraint(m, e19, 0.04*x[2]-b[23] <= 0.0)
    @constraint(m, e20, 0.0222222222222222*x[3]-b[24] <= 0.0)
    @constraint(m, e21, 0.0025*x[7]-b[22] <= 0.0)
    @constraint(m, e22, 0.0025*x[8]-b[23] <= 0.0)
    @constraint(m, e23, 0.0025*x[9]-b[24] <= 0.0)
    @constraint(m, e24, 0.00285714285714286*x[10]-b[22] <= 0.0)
    @constraint(m, e25, 0.00285714285714286*x[11]-b[23] <= 0.0)
    @constraint(m, e26, 0.00285714285714286*x[12]-b[24] <= 0.0)
    @constraint(m, e27, x[13]-b[22] <= 0.0)
    @constraint(m, e28, x[14]-b[23] <= 0.0)
    @constraint(m, e29, x[15]-b[24] <= 0.0)
    @constraint(m, e30, i[16]-3*b[22] <= 0.0)
    @constraint(m, e31, i[17]-3*b[23] <= 0.0)
    @constraint(m, e32, i[18]-3*b[24] <= 0.0)
    @constraint(m, e33, i[19]-3*b[22] <= 0.0)
    @constraint(m, e34, i[20]-3*b[23] <= 0.0)
    @constraint(m, e35, i[21]-3*b[24] <= 0.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function ex1264a(;solver = nothing)
    # Source: http://minlplib.org/ex1264a.html
    # GOpt: 8.6

    m = Model(solver)

    @variable(m, objvar)
    b_Idx = Any[17, 18, 19, 20]
    @variable(m, 0 <= b[b_Idx] <= 1, Bin)
    i_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24]
    @variable(m, i[i_Idx] >= 0, Int)

    for j=1:16
        JuMP.set_upper_bound(i[j], 5.0)
    end    
    JuMP.set_upper_bound(i[21], 15.0)
    JuMP.set_upper_bound(i[22], 12.0)
    JuMP.set_upper_bound(i[23], 9.0)
    JuMP.set_upper_bound(i[24], 6.0)

    @constraint(m, e1, -0.1*b[17]-0.2*b[18]-0.3*b[19]-0.4*b[20]-i[21]-i[22]-i[23]-i[24]+objvar == 0.0)
    @NLconstraint(m, e2, i[21]*i[1]+i[22]*i[2]+i[23]*i[3]+i[24]*i[4] >= 9.0)
    @NLconstraint(m, e3, i[21]*i[5]+i[22]*i[6]+i[23]*i[7]+i[24]*i[8] >= 7.0)
    @NLconstraint(m, e4, i[21]*i[9]+i[22]*i[10]+i[23]*i[11]+i[24]*i[12] >= 12.0)
    @NLconstraint(m, e5, i[21]*i[13]+i[22]*i[14]+i[23]*i[15]+i[24]*i[16] >= 11.0)
    @constraint(m, e6, -330*i[1]-360*i[5]-385*i[9]-415*i[13]+1700*b[17] <= 0.0)
    @constraint(m, e7, -330*i[2]-360*i[6]-385*i[10]-415*i[14]+1700*b[18] <= 0.0)
    @constraint(m, e8, -330*i[3]-360*i[7]-385*i[11]-415*i[15]+1700*b[19] <= 0.0)
    @constraint(m, e9, -330*i[4]-360*i[8]-385*i[12]-415*i[16]+1700*b[20] <= 0.0)
    @constraint(m, e10, 330*i[1]+360*i[5]+385*i[9]+415*i[13]-1900*b[17] <= 0.0)
    @constraint(m, e11, 330*i[2]+360*i[6]+385*i[10]+415*i[14]-1900*b[18] <= 0.0)
    @constraint(m, e12, 330*i[3]+360*i[7]+385*i[11]+415*i[15]-1900*b[19] <= 0.0)
    @constraint(m, e13, 330*i[4]+360*i[8]+385*i[12]+415*i[16]-1900*b[20] <= 0.0)
    @constraint(m, e14, -i[1]-i[5]-i[9]-i[13]+b[17] <= 0.0)
    @constraint(m, e15, -i[2]-i[6]-i[10]-i[14]+b[18] <= 0.0)
    @constraint(m, e16, -i[3]-i[7]-i[11]-i[15]+b[19] <= 0.0)
    @constraint(m, e17, -i[4]-i[8]-i[12]-i[16]+b[20] <= 0.0)
    @constraint(m, e18, i[1]+i[5]+i[9]+i[13]-5*b[17] <= 0.0)
    @constraint(m, e19, i[2]+i[6]+i[10]+i[14]-5*b[18] <= 0.0)
    @constraint(m, e20, i[3]+i[7]+i[11]+i[15]-5*b[19] <= 0.0)
    @constraint(m, e21, i[4]+i[8]+i[12]+i[16]-5*b[20] <= 0.0)
    @constraint(m, e22, b[17]-i[21] <= 0.0)
    @constraint(m, e23, b[18]-i[22] <= 0.0)
    @constraint(m, e24, b[19]-i[23] <= 0.0)
    @constraint(m, e25, b[20]-i[24] <= 0.0)
    @constraint(m, e26, -15*b[17]+i[21] <= 0.0)
    @constraint(m, e27, -12*b[18]+i[22] <= 0.0)
    @constraint(m, e28, -9*b[19]+i[23] <= 0.0)
    @constraint(m, e29, -6*b[20]+i[24] <= 0.0)
    @constraint(m, e30, i[21]+i[22]+i[23]+i[24] >= 8.0)
    @constraint(m, e31, -b[17]+b[18] <= 0.0)
    @constraint(m, e32, -b[18]+b[19] <= 0.0)
    @constraint(m, e33, -b[19]+b[20] <= 0.0)
    @constraint(m, e34, -i[21]+i[22] <= 0.0)
    @constraint(m, e35, -i[22]+i[23] <= 0.0)
    @constraint(m, e36, -i[23]+i[24] <= 0.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function prob03(;solver = nothing)
    # Source: http://minlplib.org/prob03.html
    # GOpt: 10.0000

    m = Model(solver)

    @variable(m, objvar)
    i_Idx = Any[1, 2]
    @variable(m, 1 <= i[i_Idx] <= 5, Int)

    @constraint(m, e1, -3*i[1] - 2*i[2] + objvar == 0.0)
    @NLconstraint(m, e2, -i[1]*i[2] <= -3.5)

    @objective(m, Min, objvar)

    return m
end

function prob10(;solver = nothing)
    # Source: http://minlplib.org/prob10.html
    # GOpt: 3.44550379

    m = Model(solver)

    @variable(m, objvar)
    x_Idx = Any[2]
    @variable(m, 0 <= x[x_Idx] <= 10)
    i_Idx = Any[3]
    @variable(m, 0 <= i[i_Idx] <= 10, Int)

    @constraint(m, e1, 0.7*x[2]+i[3] <= 7.0)
    @constraint(m, e2, 2.5*x[2]+i[3] <= 19.0)
    @NLconstraint(m, e3, 1.1*( (2*x[2]-10)^2+ (i[3]-5)^2) + sin((2*x[2]-10)^2+ (i[3]-5)^2) - objvar == 0.0)

    @objective(m, Min, objvar)

    return m
end

function st_miqp1(;solver = nothing)
    # source: http://minlplib.org/st_miqp1.html
    # GOpt: 281.0000

    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    i_Idx = Any[1, 2, 3, 4, 5]

    @variable(m, 0 <= i[i_Idx] <= 1, Int)

    # ----- Constraints ----- #
    @constraint(m, e1, 20*i[1] + 12*i[2] + 11*i[3] + 7*i[4] + 4*i[5] >= 40.0)
    @NLconstraint(m, e2, -(50*i[1]^2 + 42*i[1] + 50*i[2]^2 + 44*i[2]+ 50*i[3]^2 + 45*i[3] + 50*i[4]^2 + 47*i[4]+ 50*i[5]^2 + 47.5*i[5]) + objvar == 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function st_miqp2(;solver = nothing)
    # Source: http://minlplib.org/st_miqp2.html
    # GOpt: 2.0000
    
    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    i_Idx = Any[1, 2, 3, 4]
    @variable(m, i[i_Idx] >= 0, Int)
    JuMP.set_upper_bound(i[1], 1.0)
    JuMP.set_upper_bound(i[2], 1.0)
    JuMP.set_upper_bound(i[3], 1.0e10)
    JuMP.set_upper_bound(i[4], 1.0e10)

    # ----- Constraints ----- #
    @constraint(m, e1, -10*i[1]+i[3] <= 0.0)
    @constraint(m, e2, -20*i[2]+i[4] <= 0.0)
    @constraint(m, e3, i[3]+i[4] >= 5.0)
    @NLconstraint(m, e4, -(4*i[3]^2 - 3*i[3] + 2*i[4]^2 - 10*i[4]) - 4*i[1] -5*i[2] + objvar == 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function st_miqp3(;solver = nothing)
    # Source: http://minlplib.org/st_miqp3.html
    # GOpt: -6.0000

    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    i_Idx = Any[1, 2]
    @variable(m, i[i_Idx] >= 0, Int)
    JuMP.set_upper_bound(i[1], 3.0)
    JuMP.set_upper_bound(i[2], 1.0e15)

    # ----- Constraints ----- #
    @constraint(m, e1, -4*i[1]+i[2] <= 0.0)
    @NLconstraint(m, e2, -6*i[1]^2 + 3*i[2] + objvar == 0.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function st_miqp4(;solver = nothing)
    # Source: http://minlplib.org/st_miqp4.html
    # GOpt: -4574.0000

    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[4, 5, 6]
    @variable(m, 0 <= x[x_Idx] <= 1E15)
    i_Idx = Any[1, 2, 3]
    @variable(m, 0 <= i[i_Idx] <= 1, Int)

    # ----- Constraints ----- #
    @constraint(m, e1, x[4]+x[5]-x[6] >= 0.0)
    @constraint(m, e2, -5*i[1]+x[4] <= 0.0)
    @constraint(m, e3, -10*i[2]+x[5] <= 0.0)
    @constraint(m, e4, -30*i[3]+x[6] <= 0.0)
    @NLconstraint(m, e5, -(5*x[4]^2 + 2*x[4] + 5*x[5]^2 + 3*x[5] + 10*x[6]^2 - 500*x[6]) - 10*i[1] + 4*i[2] - 5*i[3] + objvar == 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function st_miqp5(;solver = nothing)
    # Source: http://minlplib.org/st_miqp5.html
    # GOpt: -333.88888890

    m = Model(solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[3, 4, 5, 6, 7]
    @variable(m, x[x_Idx])
    i_Idx = Any[1, 2]
    @variable(m, 0 <= i[i_Idx] <= 1, Int)
    
    JuMP.set_lower_bound(x[3], -7.24380468458)
    JuMP.set_upper_bound(x[3], 22.6826188429)
    JuMP.set_lower_bound(x[4], -6.0023781122)
    JuMP.set_upper_bound(x[4], 3.80464419615)
    JuMP.set_lower_bound(x[5], -0.797166188733)
    JuMP.set_upper_bound(x[5], 11.5189336042)
    JuMP.set_lower_bound(x[6], -8.75189948987)
    JuMP.set_upper_bound(x[6], 14.5864991498)
    JuMP.set_lower_bound(x[7], 0)
    JuMP.set_upper_bound(x[7], 19.4187214575)

    # ----- Constraints ----- #
    @constraint(m, e1, -1.93414531698*x[3]+1.80314509442*x[4]+2.89695789508*x[5]+0.729324957489*x[6]+3.8837442915*x[7] <= 60.0)
    @constraint(m, e2, -1.13150591228*x[3]+1.10500971967*x[4]-1.01838569726*x[5]+2.62556984696*x[6]+4.85468036438*x[7] <= 60.0)
    @constraint(m, e3, -0.0524800119769*x[3]-0.904837825133*x[4]+0.209520819817*x[5]-0.291729982996*x[6]-0.222506183367*x[7] <= 0.0)
    @constraint(m, e4, 0.0524800119769*x[3]+0.904837825133*x[4]-0.209520819817*x[5]+0.291729982996*x[6]+0.222506183367*x[7] <= 1.0)
    @constraint(m, e5, 0.445391966818*x[3]+0.301519984248*x[4]+0.587645368916*x[5]-0.145864991498*x[6]-0.586607210695*x[7] <= 0.0)
    @constraint(m, e6, -0.445391966818*x[3]-0.301519984248*x[4]-0.587645368916*x[5]+0.145864991498*x[6]+0.586607210695*x[7] <= 1.0)
    @constraint(m, e7, -0.328188665272*x[3]+0.199986646277*x[4]+0.506106406938*x[5]-0.583459965992*x[6]+0.505695871289*x[7] >= 0.0)
    @constraint(m, e8, -0.345682002598*x[3]-0.101625962101*x[4]+0.57594668021*x[5]+0.729324957489*x[6]+0.0809113394063*x[7] >= 0.0)
    @constraint(m, e9, 0.756087294764*x[3]-0.200079270407*x[4]+0.151379235251*x[5]+0.145864991498*x[6]+0.586607210695*x[7] >= 0.0)
    @constraint(m, e10, -i[1]+0.0524800119769*x[3]+0.904837825133*x[4]-0.209520819817*x[5]+0.291729982996*x[6]+0.222506183367*x[7] <= 0.0)
    @constraint(m, e11, i[1]-0.0524800119769*x[3]-0.904837825133*x[4]+0.209520819817*x[5]-0.291729982996*x[6]-0.222506183367*x[7] <= 0.0)
    @constraint(m, e12, -i[2]-0.445391966818*x[3]-0.301519984248*x[4]-0.587645368916*x[5]+0.145864991498*x[6]+0.586607210695*x[7] <= 0.0)
    @constraint(m, e13, i[2]+0.445391966818*x[3]+0.301519984248*x[4]+0.587645368916*x[5]-0.145864991498*x[6]-0.586607210695*x[7] <= 0.0)
    @NLconstraint(m, e14, -(5*x[6]^2 - 0.875189948987*x[6] + 52*x[7]^2 - 192.710582631*x[7]) + 54.0615511462*x[3] + 45.2691026456*x[4] + 33.0896119339*x[5] + objvar == 0.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end
