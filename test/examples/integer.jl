function ex1225a(;solver=nothing)

    m = Model(solver=solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    @variable(m, x[x_Idx])
    b_Idx = Any[22, 23, 24]
    @variable(m, b[b_Idx])
    i_Idx = Any[16, 17, 18, 19, 20, 21]
    @variable(m, i[i_Idx])
    setcategory(i[21], :Int)
    setlowerbound(i[21], 0.0)
    setupperbound(i[21], 100.0)
    setlowerbound(x[4], 0.0)
    setlowerbound(x[6], 0.0)
    setlowerbound(x[14], 0.0)
    setlowerbound(x[3], 0.0)
    setcategory(i[19], :Int)
    setlowerbound(i[19], 0.0)
    setupperbound(i[19], 100.0)
    setlowerbound(x[11], 0.0)
    setlowerbound(x[12], 0.0)
    setlowerbound(x[5], 0.0)
    setlowerbound(x[2], 0.0)
    setcategory(b[24], :Bin)
    setcategory(i[16], :Int)
    setlowerbound(i[16], 0.0)
    setupperbound(i[16], 100.0)
    setcategory(b[23], :Bin)
    setcategory(i[17], :Int)
    setlowerbound(i[17], 0.0)
    setupperbound(i[17], 100.0)
    setlowerbound(x[9], 0.0)
    setlowerbound(x[15], 0.0)
    setlowerbound(x[1], 0.0)
    setlowerbound(x[7], 0.0)
    setlowerbound(x[8], 0.0)
    setlowerbound(x[13], 0.0)
    setcategory(i[20], :Int)
    setlowerbound(i[20], 0.0)
    setupperbound(i[20], 100.0)
    setlowerbound(x[10], 0.0)
    setcategory(i[18], :Int)
    setlowerbound(i[18], 0.0)
    setupperbound(i[18], 100.0)
    setcategory(b[22], :Bin)
    setupperbound(x[1], 80.0)
    setupperbound(x[2], 25.0)
    setupperbound(x[3], 45.0)
    setupperbound(x[4], 2950.0)
    setupperbound(x[5], 2950.0)
    setupperbound(x[6], 2950.0)
    setupperbound(x[7], 400.0)
    setupperbound(x[8], 400.0)
    setupperbound(x[9], 400.0)
    setupperbound(x[10], 350.0)
    setupperbound(x[11], 350.0)
    setupperbound(x[12], 350.0)
    setupperbound(x[13], 1.0)
    setupperbound(x[14], 1.0)
    setupperbound(x[15], 1.0)
    setupperbound(i[16], 3.0)
    setupperbound(i[17], 3.0)
    setupperbound(i[18], 3.0)
    setupperbound(i[19], 3.0)
    setupperbound(i[20], 3.0)
    setupperbound(i[21], 3.0)


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

function prob10(;solver=nothing)

    m = Model(solver=solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[2]
    @variable(m, x[x_Idx])
    i_Idx = Any[3]
    @variable(m, i[i_Idx])
    setcategory(i[3], :Int)
    setlowerbound(i[3], 0.0)
    setupperbound(i[3], 100.0)
    setlowerbound(x[2], 0.0)
    setupperbound(x[2], 10.0)
    setupperbound(i[3], 10.0)


    # ----- Constraints ----- #
    @constraint(m, e1, 0.7*x[2]+i[3] <= 7.0)
    @constraint(m, e2, 2.5*x[2]+i[3] <= 19.0)
    @NLconstraint(m, e3, 1.1*( (2*x[2]-10)^2+ (i[3]-5)^2)+sin( (2*x[2]-10)^2+ (i[3]-5)^2)-objvar == 0.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function prob03(;solver=nothing)

    m = Model(solver=solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    i_Idx = Any[1, 2]
    @variable(m, i[i_Idx])
    setcategory(i[1], :Int)
    setlowerbound(i[1], 0.0)
    setupperbound(i[1], 100.0)
    setcategory(i[2], :Int)
    setlowerbound(i[2], 0.0)
    setupperbound(i[2], 100.0)
    setlowerbound(i[1], 1.0)
    setupperbound(i[1], 5.0)
    setlowerbound(i[2], 1.0)
    setupperbound(i[2], 5.0)


    # ----- Constraints ----- #
    @constraint(m, e1, -3*i[1]-2*i[2]+objvar == 0.0)
    @NLconstraint(m, e2, -i[1]*i[2] <= -3.5)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end

function st_miqp5(;solver=nothing)

    m = Model(solver=solver)

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[3, 4, 5, 6, 7]
    @variable(m, x[x_Idx])
    i_Idx = Any[1, 2]
    @variable(m, i[i_Idx])
    setcategory(i[1], :Int)
    setlowerbound(i[1], 0.0)
    setupperbound(i[1], 100.0)
    setcategory(i[2], :Int)
    setlowerbound(i[2], 0.0)
    setupperbound(i[2], 100.0)
    setupperbound(i[1], 1.0)
    setupperbound(i[2], 1.0)
    setlowerbound(x[3], -7.24380468458)
    setupperbound(x[3], 22.6826188429)
    setlowerbound(x[4], -6.0023781122)
    setupperbound(x[4], 3.80464419615)
    setlowerbound(x[5], -0.797166188733)
    setupperbound(x[5], 11.5189336042)
    setlowerbound(x[6], -8.75189948987)
    setupperbound(x[6], 14.5864991498)
    setlowerbound(x[7], 8.98296319621e-17)
    setupperbound(x[7], 19.4187214575)


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
    @NLconstraint(m, e14, -(5*x[6]*x[6]-0.875189948987*x[6]+52*x[7]*x[7]-192.710582631*x[7])+54.0615511462*x[3]+45.2691026456*x[4]+33.0896119339*x[5]+objvar == 0.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    return m
end
