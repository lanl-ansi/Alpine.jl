using POD, JuMP, Gurobi, Ipopt, MathProgBase, AmplNLWriter, CoinOptServices

function ex1264(;verbose=false, solver=nothing)

    if solver == nothing
        m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(["bonmin.iteration_limit=100"]),
                                    mip_solver=GurobiSolver(OutputFlag=0),
                                    discretization_ratio=32,
                                    log_level=100,
                                    rel_gap=0.001))
    else
        m = Model(solver=solver)
    end

    @variable(m, x[1:88])
    for i=17:68
      setcategory(x[i], :Bin)
    end
    for i=73:88
      setcategory(x[i], :Bin)
    end

    for i=1:16
      @constraint(m, x[i] >= 0.0)
      @constraint(m, x[i] <= 5.0)
    end
    @constraint(m, x[69] >= 0.0)
    @constraint(m, x[69] <= 15.0)
    @constraint(m, x[70] >= 0.0)
    @constraint(m, x[70] <= 12.0)
    @constraint(m, x[71] >= 0.0)
    @constraint(m, x[71] <= 9.0)
    @constraint(m, x[72] >= 0.0)
    @constraint(m, x[72] <= 6.0)

    @objective(m, Min, 0.1*x[65] + 0.2*x[66] + 0.3*x[67] + 0.4*x[68] + x[69] + x[70] + x[71] + x[72])
    @NLconstraint(m, x[69]*x[1]+x[70]*x[2]+x[71]*x[3]+x[72]*x[4]>=9)  #= e2: =#
    @NLconstraint(m, x[69]*x[5]+x[70]*x[6]+x[71]*x[7]+x[72]*x[8]>=7)  #= e3: =#
    @NLconstraint(m, x[69]*x[9]+x[70]*x[10]+x[71]*x[11]+x[72]*x[12]>=12)  #= e4: =#
    @NLconstraint(m, x[69]*x[13]+x[70]*x[14]+x[71]*x[15]+x[72]*x[16]>=11)  #= e5: =#

    @constraint(m, -330*x[1]-360*x[5]-385*x[9]-415*x[13]+1700*x[65]<=0)  #= e6: =#
    @constraint(m, -330*x[2]-360*x[6]-385*x[10]-415*x[14]+1700*x[66]<=0)  #= e7: =#
    @constraint(m, -330*x[3]-360*x[7]-385*x[11]-415*x[15]+1700*x[67]<=0)  #= e8: =#
    @constraint(m, -330*x[4]-360*x[8]-385*x[12]-415*x[16]+1700*x[68]<=0)  #= e9: =#
    @constraint(m, 330*x[1]+360*x[5]+385*x[9]+415*x[13]-1900*x[65]<=0)  #= e10: =#
    @constraint(m, 330*x[2]+360*x[6]+385*x[10]+415*x[14]-1900*x[66]<=0)  #= e11: =#
    @constraint(m, 330*x[3]+360*x[7]+385*x[11]+415*x[15]-1900*x[67]<=0)  #= e12: =#
    @constraint(m, 330*x[4]+360*x[8]+385*x[12]+415*x[16]-1900*x[68]<=0)  #= e13: =#
    @constraint(m, -x[1]-x[5]-x[9]-x[13]+x[65]<=0)  #= e14: =#
    @constraint(m, -x[2]-x[6]-x[10]-x[14]+x[66]<=0)  #= e15: =#
    @constraint(m, -x[3]-x[7]-x[11]-x[15]+x[67]<=0)  #= e16: =#
    @constraint(m, -x[4]-x[8]-x[12]-x[16]+x[68]<=0)  #= e17: =#
    @constraint(m, x[1]+x[5]+x[9]+x[13]-5*x[65]<=0)  #= e18: =#
    @constraint(m, x[2]+x[6]+x[10]+x[14]-5*x[66]<=0)  #= e19: =#
    @constraint(m, x[3]+x[7]+x[11]+x[15]-5*x[67]<=0)  #= e20: =#
    @constraint(m, x[4]+x[8]+x[12]+x[16]-5*x[68]<=0)  #= e21: =#
    @constraint(m, x[65]-x[69]<=0)  #= e22: =#
    @constraint(m, x[66]-x[70]<=0)  #= e23: =#
    @constraint(m, x[67]-x[71]<=0)  #= e24: =#
    @constraint(m, x[68]-x[72]<=0)  #= e25: =#
    @constraint(m, -15*x[65]+x[69]<=0)  #= e26: =#
    @constraint(m, -12*x[66]+x[70]<=0)  #= e27: =#
    @constraint(m, -9*x[67]+x[71]<=0)  #= e28: =#
    @constraint(m, -6*x[68]+x[72]<=0)  #= e29: =#
    @constraint(m, x[69]+x[70]+x[71]+x[72]>=8)  #= e30: =#
    @constraint(m, -x[65]+x[66]<=0)  #= e31: =#
    @constraint(m, -x[66]+x[67]<=0)  #= e32: =#
    @constraint(m, -x[67]+x[68]<=0)  #= e33: =#
    @constraint(m, -x[69]+x[70]<=0)  #= e34: =#
    @constraint(m, -x[70]+x[71]<=0)  #= e35: =#
    @constraint(m, -x[71]+x[72]<=0)  #= e36: =#
    @constraint(m, x[1]-x[17]-2*x[18]-4*x[19] ==0)  #= e37: =#
    @constraint(m, x[2]-x[20]-2*x[21]-4*x[22] ==0)  #= e38: =#
    @constraint(m, x[3]-x[23]-2*x[24]-4*x[25] ==0)  #= e39: =#
    @constraint(m, x[4]-x[26]-2*x[27]-4*x[28] ==0)  #= e40: =#
    @constraint(m, x[5]-x[29]-2*x[30]-4*x[31] ==0)  #= e41: =#
    @constraint(m, x[6]-x[32]-2*x[33]-4*x[34] ==0)  #= e42: =#
    @constraint(m, x[7]-x[35]-2*x[36]-4*x[37] ==0)  #= e43: =#
    @constraint(m, x[8]-x[38]-2*x[39]-4*x[40] ==0)  #= e44: =#
    @constraint(m, x[9]-x[41]-2*x[42]-4*x[43] ==0)  #= e45: =#
    @constraint(m, x[10]-x[44]-2*x[45]-4*x[46] ==0)  #= e46: =#
    @constraint(m, x[11]-x[47]-2*x[48]-4*x[49] ==0)  #= e47: =#
    @constraint(m, x[12]-x[50]-2*x[51]-4*x[52] ==0)  #= e48: =#
    @constraint(m, x[13]-x[53]-2*x[54]-4*x[55] ==0)  #= e49: =#
    @constraint(m, x[14]-x[56]-2*x[57]-4*x[58] ==0)  #= e50: =#
    @constraint(m, x[15]-x[59]-2*x[60]-4*x[61] ==0)  #= e51: =#
    @constraint(m, x[16]-x[62]-2*x[63]-4*x[64] ==0)  #= e52: =#
    @constraint(m, x[69]-x[73]-2*x[74]-4*x[75]-8*x[76] ==0)  #= e53: =#
    @constraint(m, x[70]-x[77]-2*x[78]-4*x[79]-8*x[80] ==0)  #= e54: =#
    @constraint(m, x[71]-x[81]-2*x[82]-4*x[83]-8*x[84] ==0)  #= e55: =#
    @constraint(m, x[72]-x[85]-2*x[86]-4*x[87]-8*x[88] ==0)  #= e56: =#


    if verbose
        print(m)
    end

    return m
end

function ex1265(;verbose=false, solver=nothing)

    if solver == nothing
        m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(["bonmin.iteration_limit=100"]),
                                    mip_solver=GurobiSolver(OutputFlag=0),
                                    discretization_ratio=32,
                                    log_level=100,
                                    rel_gap=0.001))
    else
        m = Model(solver=solver)
    end

    #define variables and bounds
    @variable(m, x[1:130])

    for i=26:105
      setcategory(x[i], :Bin)
    end
    for i=111:130
      setcategory(x[i], :Bin)
    end

    # Bounds Constraints
    for i=1:25
      @constraint(m, x[i] >= 0.0)
      @constraint(m, x[i] <= 5.0)
    end
    @constraint(m, x[106] >= 0.0)
    @constraint(m, x[106] <= 15.0)
    @constraint(m, x[107] >= 0.0)
    @constraint(m, x[107] <= 12.0)
    @constraint(m, x[108] >= 0.0)
    @constraint(m, x[108] <= 9.0)
    @constraint(m, x[109] >= 0.0)
    @constraint(m, x[109] <= 6.0)
    @constraint(m, x[110] >= 0.0)
    @constraint(m, x[110] <= 6.0)

    @objective(m, Min, 0.1*x[101] + 0.2*x[102] + 0.3*x[103] + 0.4*x[104] + 0.5*x[105] + x[106] + x[107] + x[108] + x[109] + x[110])
    @NLconstraint(m, x[106]*x[1]+x[107]*x[2]+x[108]*x[3]+x[109]*x[4]+x[110]*x[5]>=12)  #= e2: =#
    @NLconstraint(m, x[106]*x[6]+x[107]*x[7]+x[108]*x[8]+x[109]*x[9]+x[110]*x[10]>=6)  #= e3: =#
    @NLconstraint(m, x[106]*x[11]+x[107]*x[12]+x[108]*x[13]+x[109]*x[14]+x[110]*x[15]>=15)  #= e4: =#
    @NLconstraint(m, x[106]*x[16]+x[107]*x[17]+x[108]*x[18]+x[109]*x[19]+x[110]*x[20]>=6)  #= e5: =#
    @NLconstraint(m, x[106]*x[21]+x[107]*x[22]+x[108]*x[23]+x[109]*x[24]+x[110]*x[25]>=9)  #= e6: =#

    @constraint(m, -330*x[1]-360*x[6]-370*x[11]-415*x[16]-435*x[21]+1800*x[101]<=0)  #= e7: =#
    @constraint(m, -330*x[2]-360*x[7]-370*x[12]-415*x[17]-435*x[22]+1800*x[102]<=0)  #= e8: =#
    @constraint(m, -330*x[3]-360*x[8]-370*x[13]-415*x[18]-435*x[23]+1800*x[103]<=0)  #= e9: =#
    @constraint(m, -330*x[4]-360*x[9]-370*x[14]-415*x[19]-435*x[24]+1800*x[104]<=0)  #= e10: =#
    @constraint(m, -330*x[5]-360*x[10]-370*x[15]-415*x[20]-435*x[25]+1800*x[105]<=0)  #= e11: =#
    @constraint(m, 330*x[1]+360*x[6]+370*x[11]+415*x[16]+435*x[21]-2000*x[101]<=0)  #= e12: =#
    @constraint(m, 330*x[2]+360*x[7]+370*x[12]+415*x[17]+435*x[22]-2000*x[102]<=0)  #= e13: =#
    @constraint(m, 330*x[3]+360*x[8]+370*x[13]+415*x[18]+435*x[23]-2000*x[103]<=0)  #= e14: =#
    @constraint(m, 330*x[4]+360*x[9]+370*x[14]+415*x[19]+435*x[24]-2000*x[104]<=0)  #= e15: =#
    @constraint(m, 330*x[5]+360*x[10]+370*x[15]+415*x[20]+435*x[25]-2000*x[105]<=0)  #= e16: =#
    @constraint(m, -x[1]-x[6]-x[11]-x[16]-x[21]+x[101]<=0)  #= e17: =#
    @constraint(m, -x[2]-x[7]-x[12]-x[17]-x[22]+x[102]<=0)  #= e18: =#
    @constraint(m, -x[3]-x[8]-x[13]-x[18]-x[23]+x[103]<=0)  #= e19: =#
    @constraint(m, -x[4]-x[9]-x[14]-x[19]-x[24]+x[104]<=0)  #= e20: =#
    @constraint(m, -x[5]-x[10]-x[15]-x[20]-x[25]+x[105]<=0)  #= e21: =#
    @constraint(m, x[1]+x[6]+x[11]+x[16]+x[21]-5*x[101]<=0)  #= e22: =#
    @constraint(m, x[2]+x[7]+x[12]+x[17]+x[22]-5*x[102]<=0)  #= e23: =#
    @constraint(m, x[3]+x[8]+x[13]+x[18]+x[23]-5*x[103]<=0)  #= e24: =#
    @constraint(m, x[4]+x[9]+x[14]+x[19]+x[24]-5*x[104]<=0)  #= e25: =#
    @constraint(m, x[5]+x[10]+x[15]+x[20]+x[25]-5*x[105]<=0)  #= e26: =#
    @constraint(m, x[101]-x[106]<=0)  #= e27: =#
    @constraint(m, x[102]-x[107]<=0)  #= e28: =#
    @constraint(m, x[103]-x[108]<=0)  #= e29: =#
    @constraint(m, x[104]-x[109]<=0)  #= e30: =#
    @constraint(m, x[105]-x[110]<=0)  #= e31: =#
    @constraint(m, -15*x[101]+x[106]<=0)  #= e32: =#
    @constraint(m, -12*x[102]+x[107]<=0)  #= e33: =#
    @constraint(m, -9*x[103]+x[108]<=0)  #= e34: =#
    @constraint(m, -6*x[104]+x[109]<=0)  #= e35: =#
    @constraint(m, -6*x[105]+x[110]<=0)  #= e36: =#
    @constraint(m, x[106]+x[107]+x[108]+x[109]+x[110]>=10)  #= e37: =#
    @constraint(m, -x[101]+x[102]<=0)  #= e38: =#
    @constraint(m, -x[102]+x[103]<=0)  #= e39: =#
    @constraint(m, -x[103]+x[104]<=0)  #= e40: =#
    @constraint(m, -x[104]+x[105]<=0)  #= e41: =#
    @constraint(m, -x[106]+x[107]<=0)  #= e42: =#
    @constraint(m, -x[107]+x[108]<=0)  #= e43: =#
    @constraint(m, -x[108]+x[109]<=0)  #= e44: =#
    @constraint(m, -x[109]+x[110]<=0)  #= e45: =#
    @constraint(m, x[1]-x[26]-2*x[27]-4*x[28] ==0)  #= e46: =#
    @constraint(m, x[2]-x[29]-2*x[30]-4*x[31] ==0)  #= e47: =#
    @constraint(m, x[3]-x[32]-2*x[33]-4*x[34] ==0)  #= e48: =#
    @constraint(m, x[4]-x[35]-2*x[36]-4*x[37] ==0)  #= e49: =#
    @constraint(m, x[5]-x[38]-2*x[39]-4*x[40] ==0)  #= e50: =#
    @constraint(m, x[6]-x[41]-2*x[42]-4*x[43] ==0)  #= e51: =#
    @constraint(m, x[7]-x[44]-2*x[45]-4*x[46] ==0)  #= e52: =#
    @constraint(m, x[8]-x[47]-2*x[48]-4*x[49] ==0)  #= e53: =#
    @constraint(m, x[9]-x[50]-2*x[51]-4*x[52] ==0)  #= e54: =#
    @constraint(m, x[10]-x[53]-2*x[54]-4*x[55] ==0)  #= e55: =#
    @constraint(m, x[11]-x[56]-2*x[57]-4*x[58] ==0)  #= e56: =#
    @constraint(m, x[12]-x[59]-2*x[60]-4*x[61] ==0)  #= e57: =#
    @constraint(m, x[13]-x[62]-2*x[63]-4*x[64] ==0)  #= e58: =#
    @constraint(m, x[14]-x[65]-2*x[66]-4*x[67] ==0)  #= e59: =#
    @constraint(m, x[15]-x[68]-2*x[69]-4*x[70] ==0)  #= e60: =#
    @constraint(m, x[16]-x[71]-2*x[72]-4*x[73] ==0)  #= e61: =#
    @constraint(m, x[17]-x[74]-2*x[75]-4*x[76] ==0)  #= e62: =#
    @constraint(m, x[18]-x[77]-2*x[78]-4*x[79] ==0)  #= e63: =#
    @constraint(m, x[19]-x[80]-2*x[81]-4*x[82] ==0)  #= e64: =#
    @constraint(m, x[20]-x[83]-2*x[84]-4*x[85] ==0)  #= e65: =#
    @constraint(m, x[21]-x[86]-2*x[87]-4*x[88] ==0)  #= e66: =#
    @constraint(m, x[22]-x[89]-2*x[90]-4*x[91] ==0)  #= e67: =#
    @constraint(m, x[23]-x[92]-2*x[93]-4*x[94] ==0)  #= e68: =#
    @constraint(m, x[24]-x[95]-2*x[96]-4*x[97] ==0)  #= e69: =#
    @constraint(m, x[25]-x[98]-2*x[99]-4*x[100] ==0)  #= e70: =#
    @constraint(m, x[106]-x[111]-2*x[112]-4*x[113]-8*x[114] ==0)  #= e71: =#
    @constraint(m, x[107]-x[115]-2*x[116]-4*x[117]-8*x[118] ==0)  #= e72: =#
    @constraint(m, x[108]-x[119]-2*x[120]-4*x[121]-8*x[122] ==0)  #= e73: =#
    @constraint(m, x[109]-x[123]-2*x[124]-4*x[125]-8*x[126] ==0)  #= e74: =#
    @constraint(m, x[110]-x[127]-2*x[128]-4*x[129]-8*x[130] ==0)  #= e75: =#

    if verbose
        print(m)
    end

    return m
end

function ex1266(;verbose=false, solver=nothing)

    if solver == nothing
        m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(["bonmin.iteration_limit=100"]),
                                    mip_solver=GurobiSolver(OutputFlag=0),
                                    discretization_ratio=32,
                                    log_level=100,
                                    rel_gap=0.001))
    else
        m = Model(solver=solver)
    end

    #define variables and bounds
    @variable(m, x[1:180])
    for i=37:150
      setcategory(x[i], :Bin)
    end
    for i=157:180
      setcategory(x[i], :Bin)
    end

    # Bounds Constraints
    for i=1:36
      @constraint(m, x[i] >= 0.0)
      @constraint(m, x[i] <= 5.0)
    end
    @constraint(m, x[151] >= 0.0)
    @constraint(m, x[151] <= 15.0)
    @constraint(m, x[152] >= 0.0)
    @constraint(m, x[152] <= 12.0)
    @constraint(m, x[153] >= 0.0)
    @constraint(m, x[153] <= 8.0)
    @constraint(m, x[154] >= 0.0)
    @constraint(m, x[154] <= 7.0)
    @constraint(m, x[155] >= 0.0)
    @constraint(m, x[155] <= 4.0)
    @constraint(m, x[156] >= 0.0)
    @constraint(m, x[156] <= 2.0)

    @objective(m, Min, 0.1*x[145] + 0.2*x[146] + 0.3*x[147] + 0.4*x[148] + 0.5*x[149] + x[151] + x[152] + x[153] + x[154] + x[155] + x[156])
    @NLconstraint(m, x[151]*x[1]+x[152]*x[2]+x[153]*x[3]+x[154]*x[4]+x[155]*x[5]+x[156]*x[6]>=8)  #= e2: =#
    @NLconstraint(m, x[151]*x[7]+x[152]*x[8]+x[153]*x[9]+x[154]*x[10]+x[155]*x[11]+x[156]*x[12]>=16)  #= e3: =#
    @NLconstraint(m, x[151]*x[13]+x[152]*x[14]+x[153]*x[15]+x[154]*x[16]+x[155]*x[17]+x[156]*x[18]>=12)  #= e4: =#
    @NLconstraint(m, x[151]*x[19]+x[152]*x[20]+x[153]*x[21]+x[154]*x[22]+x[155]*x[23]+x[156]*x[24]>=7)  #= e5: =#
    @NLconstraint(m, x[151]*x[25]+x[152]*x[26]+x[153]*x[27]+x[154]*x[28]+x[155]*x[29]+x[156]*x[30]>=14)  #= e6: =#
    @NLconstraint(m, x[151]*x[31]+x[152]*x[32]+x[153]*x[33]+x[154]*x[34]+x[155]*x[35]+x[156]*x[36]>=16)  #= e7: =#

    @constraint(m, -330*x[1]-360*x[7]-380*x[13]-430*x[19]-490*x[25]-530*x[31]+2100*x[145]<=0)  #= e8: =#
    @constraint(m, -330*x[2]-360*x[8]-380*x[14]-430*x[20]-490*x[26]-530*x[32]+2100*x[146]<=0)  #= e9: =#
    @constraint(m, -330*x[3]-360*x[9]-380*x[15]-430*x[21]-490*x[27]-530*x[33]+2100*x[147]<=0)  #= e10: =#
    @constraint(m, -330*x[4]-360*x[10]-380*x[16]-430*x[22]-490*x[28]-530*x[34]+2100*x[148]<=0)  #= e11: =#
    @constraint(m, -330*x[5]-360*x[11]-380*x[17]-430*x[23]-490*x[29]-530*x[35]+2100*x[149]<=0)  #= e12: =#
    @constraint(m, -330*x[6]-360*x[12]-380*x[18]-430*x[24]-490*x[30]-530*x[36]+2100*x[150]<=0)  #= e13: =#
    @constraint(m, 330*x[1]+360*x[7]+380*x[13]+430*x[19]+490*x[25]+530*x[31]-2200*x[145]<=0)  #= e14: =#
    @constraint(m, 330*x[2]+360*x[8]+380*x[14]+430*x[20]+490*x[26]+530*x[32]-2200*x[146]<=0)  #= e15: =#
    @constraint(m, 330*x[3]+360*x[9]+380*x[15]+430*x[21]+490*x[27]+530*x[33]-2200*x[147]<=0)  #= e16: =#
    @constraint(m, 330*x[4]+360*x[10]+380*x[16]+430*x[22]+490*x[28]+530*x[34]-2200*x[148]<=0)  #= e17: =#
    @constraint(m, 330*x[5]+360*x[11]+380*x[17]+430*x[23]+490*x[29]+530*x[35]-2200*x[149]<=0)  #= e18: =#
    @constraint(m, 330*x[6]+360*x[12]+380*x[18]+430*x[24]+490*x[30]+530*x[36]-2200*x[150]<=0)  #= e19: =#
    @constraint(m, -x[1]-x[7]-x[13]-x[19]-x[25]-x[31]+x[145]<=0)  #= e20: =#
    @constraint(m, -x[2]-x[8]-x[14]-x[20]-x[26]-x[32]+x[146]<=0)  #= e21: =#
    @constraint(m, -x[3]-x[9]-x[15]-x[21]-x[27]-x[33]+x[147]<=0)  #= e22: =#
    @constraint(m, -x[4]-x[10]-x[16]-x[22]-x[28]-x[34]+x[148]<=0)  #= e23: =#
    @constraint(m, -x[5]-x[11]-x[17]-x[23]-x[29]-x[35]+x[149]<=0)  #= e24: =#
    @constraint(m, -x[6]-x[12]-x[18]-x[24]-x[30]-x[36]+x[150]<=0)  #= e25: =#
    @constraint(m, x[1]+x[7]+x[13]+x[19]+x[25]+x[31]-5*x[145]<=0)  #= e26: =#
    @constraint(m, x[2]+x[8]+x[14]+x[20]+x[26]+x[32]-5*x[146]<=0)  #= e27: =#
    @constraint(m, x[3]+x[9]+x[15]+x[21]+x[27]+x[33]-5*x[147]<=0)  #= e28: =#
    @constraint(m, x[4]+x[10]+x[16]+x[22]+x[28]+x[34]-5*x[148]<=0)  #= e29: =#
    @constraint(m, x[5]+x[11]+x[17]+x[23]+x[29]+x[35]-5*x[149]<=0)  #= e30: =#
    @constraint(m, x[6]+x[12]+x[18]+x[24]+x[30]+x[36]-5*x[150]<=0)  #= e31: =#
    @constraint(m, x[145]-x[151]<=0)  #= e32: =#
    @constraint(m, x[146]-x[152]<=0)  #= e33: =#
    @constraint(m, x[147]-x[153]<=0)  #= e34: =#
    @constraint(m, x[148]-x[154]<=0)  #= e35: =#
    @constraint(m, x[149]-x[155]<=0)  #= e36: =#
    @constraint(m, x[150]-x[156]<=0)  #= e37: =#
    @constraint(m, -15*x[145]+x[151]<=0)  #= e38: =#
    @constraint(m, -8*x[147]+x[153]<=0)  #= e40: =#
    @constraint(m, -7*x[148]+x[154]<=0)  #= e41: =#
    @constraint(m, -4*x[149]+x[155]<=0)  #= e42: =#
    @constraint(m, -2*x[150]+x[156]<=0)  #= e43: =#
    @constraint(m, x[151]+x[152]+x[153]+x[154]+x[155]+x[156]>=16)  #= e44: =#
    @constraint(m, -x[145]+x[146]<=0)  #= e45: =#
    @constraint(m, -x[146]+x[147]<=0)  #= e46: =#
    @constraint(m, -x[147]+x[148]<=0)  #= e47: =#
    @constraint(m, -x[148]+x[149]<=0)  #= e48: =#
    @constraint(m, -x[149]+x[150]<=0)  #= e49: =#
    @constraint(m, -x[151]+x[152]<=0)  #= e50: =#
    @constraint(m, -x[152]+x[153]<=0)  #= e51: =#
    @constraint(m, -x[153]+x[154]<=0)  #= e52: =#
    @constraint(m, -x[154]+x[155]<=0)  #= e53: =#
    @constraint(m, -x[155]+x[156]<=0)  #= e54: =#
    @constraint(m, x[1]-x[37]-2*x[38]-4*x[39] ==0)  #= e55: =#
    @constraint(m, x[2]-x[40]-2*x[41]-4*x[42] ==0)  #= e56: =#
    @constraint(m, x[3]-x[43]-2*x[44]-4*x[45] ==0)  #= e57: =#
    @constraint(m, x[4]-x[46]-2*x[47]-4*x[48] ==0)  #= e58: =#
    @constraint(m, x[5]-x[49]-2*x[50]-4*x[51] ==0)  #= e59: =#
    @constraint(m, x[6]-x[52]-2*x[53]-4*x[54] ==0)  #= e60: =#
    @constraint(m, x[7]-x[55]-2*x[56]-4*x[57] ==0)  #= e61: =#
    @constraint(m, x[8]-x[58]-2*x[59]-4*x[60] ==0)  #= e62: =#
    @constraint(m, x[9]-x[61]-2*x[62]-4*x[63] ==0)  #= e63: =#
    @constraint(m, x[10]-x[64]-2*x[65]-4*x[66] ==0)  #= e64: =#
    @constraint(m, x[11]-x[67]-2*x[68]-4*x[69] ==0)  #= e65: =#
    @constraint(m, x[12]-x[70]-2*x[71]-4*x[72] ==0)  #= e66: =#
    @constraint(m, x[13]-x[73]-2*x[74]-4*x[75] ==0)  #= e67: =#
    @constraint(m, x[14]-x[76]-2*x[77]-4*x[78] ==0)  #= e68: =#
    @constraint(m, x[15]-x[79]-2*x[80]-4*x[81] ==0)  #= e69: =#
    @constraint(m, x[16]-x[82]-2*x[83]-4*x[84] ==0)  #= e70: =#
    @constraint(m, x[17]-x[85]-2*x[86]-4*x[87] ==0)  #= e71: =#
    @constraint(m, x[18]-x[88]-2*x[89]-4*x[90] ==0)  #= e72: =#
    @constraint(m, x[19]-x[91]-2*x[92]-4*x[93] ==0)  #= e73: =#
    @constraint(m, x[20]-x[94]-2*x[95]-4*x[96] ==0)  #= e74: =#
    @constraint(m, x[21]-x[97]-2*x[98]-4*x[99] ==0)  #= e75: =#
    @constraint(m, x[22]-x[100]-2*x[101]-4*x[102] ==0)  #= e76: =#
    @constraint(m, x[23]-x[103]-2*x[104]-4*x[105] ==0)  #= e77: =#
    @constraint(m, x[24]-x[106]-2*x[107]-4*x[108] ==0)  #= e78: =#
    @constraint(m, x[25]-x[109]-2*x[110]-4*x[111] ==0)  #= e79: =#
    @constraint(m, x[26]-x[112]-2*x[113]-4*x[114] ==0)  #= e80: =#
    @constraint(m, x[27]-x[115]-2*x[116]-4*x[117] ==0)  #= e81: =#
    @constraint(m, x[28]-x[118]-2*x[119]-4*x[120] ==0)  #= e82: =#
    @constraint(m, x[29]-x[121]-2*x[122]-4*x[123] ==0)  #= e83: =#
    @constraint(m, x[30]-x[124]-2*x[125]-4*x[126] ==0)  #= e84: =#
    @constraint(m, x[31]-x[127]-2*x[128]-4*x[129] ==0)  #= e85: =#
    @constraint(m, x[32]-x[130]-2*x[131]-4*x[132] ==0)  #= e86: =#
    @constraint(m, x[33]-x[133]-2*x[134]-4*x[135] ==0)  #= e87: =#
    @constraint(m, x[34]-x[136]-2*x[137]-4*x[138] ==0)  #= e88: =#
    @constraint(m, x[35]-x[139]-2*x[140]-4*x[141] ==0)  #= e89: =#
    @constraint(m, x[36]-x[142]-2*x[143]-4*x[144] ==0)  #= e90: =#
    @constraint(m, x[151]-x[157]-2*x[158]-4*x[159]-8*x[160] ==0)  #= e91: =#
    @constraint(m, x[152]-x[161]-2*x[162]-4*x[163]-8*x[164] ==0)  #= e92: =#
    @constraint(m, x[153]-x[165]-2*x[166]-4*x[167]-8*x[168] ==0)  #= e93: =#
    @constraint(m, x[154]-x[169]-2*x[170]-4*x[171]-8*x[172] ==0)  #= e94: =#
    @constraint(m, x[155]-x[173]-2*x[174]-4*x[175]-8*x[176] ==0)  #= e95: =#
    @constraint(m, x[156]-x[177]-2*x[178]-4*x[179]-8*x[180] ==0)  #= e96: =#

    if verbose
        print(m)
    end

    return m
end

function ex1223a(;verbose=false, solver=nothing)

    if solver == nothing
        m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(["bonmin.iteration_limit=100"]),
                                    mip_solver=GurobiSolver(OutputFlag=0),
                                    discretization_ratio=32,
                                    log_level=100,
                                    rel_gap=0.001))
    else
        m = Model(solver=solver)
    end


    #define variables and bounds
    @variable(m, x[1:7])
    setcategory(x[4], :Bin)
    setcategory(x[5], :Bin)
    setcategory(x[6], :Bin)
    setcategory(x[7], :Bin)

    @constraint(m, x[1] >= 0.0)
    @constraint(m, x[1] <= 10.0)
    @constraint(m, x[2] >= 0.0)
    @constraint(m, x[2] <= 10.0)
    @constraint(m, x[3] >= 0.0)
    @constraint(m, x[3] <= 10.0)

    @constraint(m, x[1]+x[2]+x[3]+x[4]+x[5]+x[6]<=5)  #= e1: =#
    @constraint(m, x[1]+x[4]<=1.2)  #= e3: =#
    @constraint(m, x[2]+x[5]<=1.8)  #= e4: =#
    @constraint(m, x[3]+x[6]<=2.5)  #= e5: =#
    @constraint(m, x[1]+x[7]<=1.2)  #= e6: =#

    # @objective(m, Min, (-1 + x[1])^2 + (-2 + x[2])^2 + (-3 + x[3])^2 - x[4] - 3*x[5] - x[6] - 0.693147180559945*x[7] + 6)
    # Expanded objective
    @objective(m, Min, x[1]^2-2*x[1]+1 + x[2]^2-4*x[2]+4 + x[3]^2-6*x[3]+9 - x[4] - 3*x[5] - x[6] - 0.693147180559945*x[7] + 6)
    @NLconstraint(m, x[1]^2+x[2]^2+x[3]^2+x[6]<=5.5)  #= e2: =#
    @NLconstraint(m, x[2]^2+x[5]<=1.64)  #= e7: =#
    @NLconstraint(m, x[3]^2+x[6]<=4.25)  #= e8: =#
    @NLconstraint(m, x[3]^2+x[5]<=4.64)  #= e9: =#

    if verbose
        print(m)
    end

    return m
end
