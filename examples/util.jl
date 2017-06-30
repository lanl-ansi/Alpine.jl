using POD, JuMP, Gurobi, AmplNLWriter, CoinOptServices, MathProgBase

function util(;verbose=false, solver=nothing)

    if solver==nothing
        m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(),
    							   mip_solver=GurobiSolver(OutputFlag=0),
                                   log_level=1))
    else
        m = Model(solver=solver)
    end

    @variable(m, x[1:145])
    for i=1:28
        setcategory(x[i], :Bin)
    end

    lb=zeros(145)
    lb[58] = 353.66
    lb[76] = 32.8
    lb[93] = 24.1
    lb[97] = 1.2
    lb[110] = 8.7

    for i=29:145
        setlowerbound(x[i], lb[i])
    end

    setupperbound(x[47], 29.2)
    setupperbound(x[48], 1.2)
    setupperbound(x[95], 21.0)
    setupperbound(x[97], 2.2)
    setupperbound(x[105], 98.72)

    # Artifical Upper bounds
    setupperbound(x[51], 2000)
    setupperbound(x[52], 2000)
    setupperbound(x[53], 2000)
    setupperbound(x[54], 2000)
    setupperbound(x[49], 2000)
    setupperbound(x[59], 2000)
    setupperbound(x[77], 2000)
    # Artifical Upper bounds

    @objective(m, Min, 87.3969*x[29] + 0.03781*x[30] + 0.03781*x[31] + 0.03781*x[32] + 0.03781*x[33] + 0.03781*x[34] + 0.03781*x[35] + 0.03781*x[36] + 0.03781*x[37] + 0.03781*x[38] + 0.03781*x[39] + 0.03781*x[40] + 0.03781*x[41] + 0.03781*x[42] + 0.03781*x[43] + 0.3*x[44] + 0.017*x[45])

    @NLconstraint(m, x[53]*x[49]+x[52]*x[51]-0.5439*x[44]-0.5439*x[46]-1.22963*x[47]-2.74289*x[48]-6.94492*x[50] ==0)  #= e1: =#
    @NLconstraint(m, x[53]*x[54]-3.22692*x[55]-9.05971*x[56]+x[57]-x[58] ==0)  #= e2: =#
    @NLconstraint(m, x[53]*x[59]+7.32917*x[60]+7.70075*x[61]-7.31039*x[62]-7.31039*x[63]-7.31039*x[64]-7.31039*x[65]-7.31039*x[66]-7.31039*x[67]-7.31039*x[68]-7.31039*x[69]-7.31039*x[70]-7.31039*x[71]-7.31039*x[72]-7.31039*x[73]-7.31039*x[74]-7.31039*x[75]-7.31039*x[76] ==0)  #= e3: =#
    @NLconstraint(m, x[53]*x[77]-6.94492*x[50]+7.31039*x[62]+7.03739*x[78]+7.03739*x[79]+6.91604*x[80]+6.91604*x[81]+6.91604*x[82]+6.91604*x[83]+6.91604*x[84]+6.91604*x[85]+6.91604*x[86]+6.91604*x[87]+6.91604*x[88]+6.91604*x[89]+6.91604*x[90]+6.91604*x[91]-6.94492*x[92]-6.94492*x[93] == 0)  #= e4: =#

    @constraint(m, -x[44]-x[46]-x[47]-x[48]+x[49]-x[50]+x[51] ==0)  #= e5: =#
    @constraint(m, x[44]-x[94]-x[95] ==0)  #= e6: =#
    @constraint(m, -0.08*x[50]+x[51] ==0)  #= e7: =#
    @constraint(m, 0.05391*x[96]-x[97] ==4.45329)  #= e8: =#
    @constraint(m, x[52]-0.34851*x[53] ==6.04388)  #= e9: =#
    @constraint(m, x[53]-0.18673*x[97] ==0.82639)  #= e10: =#
    @constraint(m, x[49]-x[54]-x[59]-x[77]-x[98] ==0)  #= e11: =#
    @constraint(m, x[98]>=104.21999)  #= e12: =#
    @constraint(m, -0.5*x[49]+x[99] ==0)  #= e13: =#
    @constraint(m, -1.89474*x[99]+x[100] ==230)  #= e14: =#
    @constraint(m, x[100]-x[101] ==0)  #= e15: =#
    @constraint(m, x[54]-x[55]-x[56] ==0)  #= e16: =#
    @constraint(m, -0.05*x[54]+x[55] ==0)  #= e17: =#
    @constraint(m, -116.41403*x[29]+x[57] ==0)  #= e18: =#
    @constraint(m, x[29]-x[102] ==4.035)  #= e19: =#
    @constraint(m, -23.39292*x[102]+x[103] ==0)  #= e20: =#
    @constraint(m, -3.61341*x[103]+x[104] ==0)  #= e21: =#
    @constraint(m, x[56]-x[61]+x[105]-x[106]-x[107]-x[108]-x[109]-x[110] ==0)  #= e22: =#
    @constraint(m, x[59]+x[60]+x[61]-x[62]-x[63]-x[64]-x[65]-x[66]-x[67]-x[68]-x[69]-x[70]-x[71]-x[72]-x[73]-x[74]-x[75]-x[76] ==0)  #= e23: =#
    @constraint(m, -x[50]+x[62]+x[77]+x[78]+x[79]+x[80]+x[81]+x[82]+x[83]+x[84]+x[85]+x[86]+x[87]+x[88]+x[89]+x[90]+x[91]-x[92]-x[93] ==0)  #= e24: =#
    @constraint(m, -x[46]+x[111]+x[112] ==0)  #= e25: =#
    @constraint(m, x[106]-x[113] ==0)  #= e26: =#
    @constraint(m, -x[60]+x[107] ==0)  #= e27: =#
    @constraint(m, -x[78]+x[108] ==0)  #= e28: =#
    @constraint(m, -x[79]+x[109] ==0)  #= e29: =#
    @constraint(m, 7.70075*x[106]-6.08959*x[113]-0.00641*x[114] ==0)  #= e30: =#
    @constraint(m, -7.32917*x[60]+7.70075*x[107]-0.00641*x[115] ==0)  #= e31: =#
    @constraint(m, -7.03739*x[78]+7.70075*x[108]-0.00641*x[116] ==0)  #= e32: =#
    @constraint(m, -7.03739*x[79]+7.70075*x[109]-0.00641*x[117] ==0)  #= e33: =#
    @constraint(m, x[63]-x[80] ==0)  #= e34: =#
    @constraint(m, x[64]-x[81] ==0)  #= e35: =#
    @constraint(m, x[65]-x[82] ==0)  #= e36: =#
    @constraint(m, x[66]-x[118] ==0)  #= e37: =#
    @constraint(m, x[67]-x[83] ==0)  #= e38: =#
    @constraint(m, x[68]-x[84] ==0)  #= e39: =#
    @constraint(m, x[69]-x[85] ==0)  #= e40: =#
    @constraint(m, x[70]-x[86] ==0)  #= e41: =#
    @constraint(m, x[71]-x[87] ==0)  #= e42: =#
    @constraint(m, x[72]-x[88] ==0)  #= e43: =#
    @constraint(m, x[73]-x[89] ==0)  #= e44: =#
    @constraint(m, x[74]-x[90] ==0)  #= e45: =#
    @constraint(m, x[75]-x[91] ==0)  #= e46: =#
    @constraint(m, 7.31039*x[63]-6.91604*x[80]-0.00641*x[119] ==0)  #= e47: =#
    @constraint(m, 7.31039*x[64]-6.91604*x[81]-0.00641*x[120] ==0)  #= e48: =#
    @constraint(m, 7.31039*x[65]-6.91604*x[82]-0.00641*x[121] ==0)  #= e49: =#
    @constraint(m, 7.31039*x[66]-5.77083*x[118]-0.00641*x[122] ==0)  #= e50: =#
    @constraint(m, 7.31039*x[67]-6.91604*x[83]-0.00641*x[123] ==0)  #= e51: =#
    @constraint(m, 7.31039*x[68]-6.91604*x[84]-0.00641*x[124] ==0)  #= e52: =#
    @constraint(m, 7.31039*x[69]-6.91604*x[85]-0.00641*x[125] ==0)  #= e53: =#
    @constraint(m, 7.31039*x[70]-6.91604*x[86]-0.00641*x[126] ==0)  #= e54: =#
    @constraint(m, 7.31039*x[71]-6.91604*x[87]-0.00641*x[127] ==0)  #= e55: =#
    @constraint(m, 7.31039*x[72]-6.91604*x[88]-0.00641*x[128] ==0)  #= e56: =#
    @constraint(m, 7.31039*x[73]-6.91604*x[89]-0.00641*x[129] ==0)  #= e57: =#
    @constraint(m, 7.31039*x[74]-6.91604*x[90]-0.00641*x[130] ==0)  #= e58: =#
    @constraint(m, 7.31039*x[75]-6.91604*x[91]-0.00641*x[131] ==0)  #= e59: =#
    @constraint(m, x[111]-x[113] ==0)  #= e60: =#
    @constraint(m, x[112]-x[118] ==0)  #= e61: =#
    @constraint(m, 5.6972*x[111]-0.12983*x[132] ==0)  #= e62: =#
    @constraint(m, 5.6972*x[112]-0.12983*x[133] ==0)  #= e63: =#
    @constraint(m, x[45]-x[132]-x[133] ==4117)  #= e64: =#
    @constraint(m, -0.33333*x[45]+x[134] ==0)  #= e65: =#
    @constraint(m, -0.20286*x[134]+x[135] ==0)  #= e66: =#
    @constraint(m, x[135]-x[136] ==0)  #= e67: =#
    @constraint(m, x[135]-x[137] ==0)  #= e68: =#
    @constraint(m, x[138] ==1)  #= e69: =#
    @constraint(m, x[139] ==1)  #= e70: =#
    @constraint(m, x[140] ==1)  #= e71: =#
    @constraint(m, x[141] ==1)  #= e72: =#
    @constraint(m, x[138]-x[142] ==0)  #= e73: =#
    @constraint(m, x[138]-x[143] ==0)  #= e74: =#
    @constraint(m, x[139]-x[144] ==0)  #= e75: =#
    @constraint(m, 0.01*x[111]-100*x[144]<=0)  #= e76: =#
    @constraint(m, x[141]-x[145] ==0)  #= e77: =#
    @constraint(m, 0.01*x[112]-100*x[145]<=0)  #= e78: =#
    @constraint(m, 0.1*x[49]-100*x[142]<=0)  #= e79: =#
    @constraint(m, 0.1*x[103]-100*x[143]<=0)  #= e80: =#
    @constraint(m, 0.1*x[56]-100*x[138]<=0)  #= e81: =#
    @constraint(m, 0.1*x[106]-100*x[139]<=0)  #= e82: =#
    @constraint(m, 0.1*x[107]-100*x[140]<=0)  #= e83: =#
    @constraint(m, -100*x[1]+0.1*x[108]<=0)  #= e84: =#
    @constraint(m, -100*x[2]+0.1*x[109]<=0)  #= e85: =#
    @constraint(m, 0.1*x[60]-100*x[140]<=0)  #= e86: =#
    @constraint(m, -100*x[1]+0.1*x[78]<=0)  #= e87: =#
    @constraint(m, -100*x[2]+0.1*x[79]<=0)  #= e88: =#
    @constraint(m, 0.1*x[113]-100*x[139]<=0)  #= e89: =#
    @constraint(m, -100*x[3]+0.1*x[63]<=0)  #= e90: =#
    @constraint(m, -100*x[4]+0.1*x[64]<=0)  #= e91: =#
    @constraint(m, -x[5]+0.1*x[65]<=0)  #= e92: =#
    @constraint(m, 0.1*x[66]-100*x[141]<=0)  #= e93: =#
    @constraint(m, -100*x[6]+0.1*x[67]<=0)  #= e94: =#
    @constraint(m, -100*x[7]+0.1*x[68]<=0)  #= e95: =#
    @constraint(m, -100*x[8]+0.1*x[69]<=0)  #= e96: =#
    @constraint(m, -100*x[9]+0.1*x[70]<=0)  #= e97: =#
    @constraint(m, -100*x[10]+0.1*x[71]<=0)  #= e98: =#
    @constraint(m, -100*x[11]+0.1*x[72]<=0)  #= e99: =#
    @constraint(m, -100*x[12]+0.1*x[73]<=0)  #= e100: =#
    @constraint(m, -100*x[13]+0.1*x[74]<=0)  #= e101: =#
    @constraint(m, -100*x[14]+0.1*x[75]<=0)  #= e102: =#
    @constraint(m, -100*x[3]+0.1*x[80]<=0)  #= e103: =#
    @constraint(m, -100*x[4]+0.1*x[81]<=0)  #= e104: =#
    @constraint(m, -100*x[5]+0.1*x[82]<=0)  #= e105: =#
    @constraint(m, -100*x[6]+0.1*x[83]<=0)  #= e106: =#
    @constraint(m, -100*x[7]+0.1*x[84]<=0)  #= e107: =#
    @constraint(m, -100*x[8]+0.1*x[85]<=0)  #= e108: =#
    @constraint(m, -100*x[9]+0.1*x[86]<=0)  #= e109: =#
    @constraint(m, -100*x[10]+0.1*x[87]<=0)  #= e110: =#
    @constraint(m, -100*x[11]+0.1*x[88]<=0)  #= e111: =#
    @constraint(m, -100*x[12]+0.1*x[89]<=0)  #= e112: =#
    @constraint(m, -100*x[13]+0.1*x[90]<=0)  #= e113: =#
    @constraint(m, -100*x[14]+0.1*x[91]<=0)  #= e114: =#
    @constraint(m, 0.1*x[118]-100*x[141]<=0)  #= e115: =#
    @constraint(m, 0.01*x[114]-163.03*x[139] ==0)  #= e116: =#
    @constraint(m, 0.01*x[115]-56*x[140] ==0)  #= e117: =#
    @constraint(m, -444*x[3]+x[119] ==0)  #= e118: =#
    @constraint(m, -100*x[4]+x[120] ==0)  #= e119: =#
    @constraint(m, -125*x[5]+x[121] ==0)  #= e120: =#
    @constraint(m, 0.01*x[122]-94.5*x[141] ==0)  #= e121: =#
    @constraint(m, -75*x[6]+x[123] ==0)  #= e122: =#
    @constraint(m, -60*x[7]+x[124] ==0)  #= e123: =#
    @constraint(m, -36*x[8]+x[125] ==0)  #= e124: =#
    @constraint(m, -26*x[9]+x[126] ==0)  #= e125: =#
    @constraint(m, -445*x[10]+x[127] ==0)  #= e126: =#
    @constraint(m, -444*x[15]+x[30] ==0)  #= e127: =#
    @constraint(m, -100*x[16]+x[31] ==0)  #= e128: =#
    @constraint(m, -125*x[17]+x[32] ==0)  #= e129: =#
    @constraint(m, -75*x[18]+x[33] ==0)  #= e130: =#
    @constraint(m, -60*x[19]+x[34] ==0)  #= e131: =#
    @constraint(m, -36*x[20]+x[35] ==0)  #= e132: =#
    @constraint(m, -26*x[21]+x[36] ==0)  #= e133: =#
    @constraint(m, -445*x[22]+x[37] ==0)  #= e134: =#
    @constraint(m, -300*x[1]+0.01*x[116]<=0)  #= e135: =#
    @constraint(m, -300*x[2]+0.01*x[117]<=0)  #= e136: =#
    @constraint(m, -300*x[11]+0.01*x[128]<=0)  #= e137: =#
    @constraint(m, -300*x[12]+0.01*x[129]<=0)  #= e138: =#
    @constraint(m, -300*x[13]+0.01*x[130]<=0)  #= e139: =#
    @constraint(m, -300*x[14]+0.01*x[131]<=0)  #= e140: =#
    @constraint(m, -300*x[23]+0.01*x[38]<=0)  #= e141: =#
    @constraint(m, -300*x[24]+0.01*x[39]<=0)  #= e142: =#
    @constraint(m, -300*x[25]+0.01*x[40]<=0)  #= e143: =#
    @constraint(m, -300*x[26]+0.01*x[41]<=0)  #= e144: =#
    @constraint(m, -300*x[27]+0.01*x[42]<=0)  #= e145: =#
    @constraint(m, -300*x[28]+0.01*x[43]<=0)  #= e146: =#
    @constraint(m, x[38]-x[100]+x[116] ==0)  #= e147: =#
    @constraint(m, x[39]-x[101]+x[117] ==0)  #= e148: =#
    @constraint(m, x[40]-x[104]+x[128] ==0)  #= e149: =#
    @constraint(m, x[41]+x[129]-x[135] ==0)  #= e150: =#
    @constraint(m, x[42]+x[130]-x[136] ==0)  #= e151: =#
    @constraint(m, x[43]+x[131]-x[137] ==0)  #= e152: =#
    @constraint(m, x[1]+x[23] ==1)  #= e153: =#
    @constraint(m, x[2]+x[24] ==1)  #= e154: =#
    @constraint(m, x[3]+x[15] ==1)  #= e155: =#
    @constraint(m, x[4]+x[16] ==1)  #= e156: =#
    @constraint(m, x[5]+x[17] ==1)  #= e157: =#
    @constraint(m, x[6]+x[18] ==1)  #= e158: =#
    @constraint(m, x[7]+x[19] ==1)  #= e159: =#
    @constraint(m, x[8]+x[20] ==1)  #= e160: =#
    @constraint(m, x[9]+x[21] ==1)  #= e161: =#
    @constraint(m, x[10]+x[22] ==1)  #= e162: =#
    @constraint(m, x[11]+x[25] ==1)  #= e163: =#
    @constraint(m, x[12]+x[26] ==1)  #= e164: =#
    @constraint(m, x[13]+x[27] ==1)  #= e165: =#
    @constraint(m, x[14]+x[28] ==1)  #= e166: =#
    @constraint(m, x[14]+x[24] ==0)  #= e167: =#

    return m
end
