function blend029(;solver=nothing)

    m = Model(solver)

    @variable(m, x[1:102])
    for i=67:102
      set_binary(x[i])
    end

    for i=1:48
      set_lower_bound(x[i], 0)
      set_upper_bound(x[i], 1)
    end
    for i=49:66
      set_lower_bound(x[i], 0)
      set_upper_bound(x[i], 2)
    end

    @objective(m, Max, - 1.74*x[1] - 1.74*x[2] - 1.74*x[3] - 1.45*x[4] - 1.45*x[5] - 1.45*x[6] + 7.38*x[7] + 7.38*x[8] + 7.38*x[9] + 5.6*x[10] + 5.6*x[11] + 5.6*x[12] - 1.7*x[13] - 1.7*x[14] - 1.7*x[15] - 1.18*x[16] - 1.18*x[17] - 1.18*x[18] + 7.21*x[19] + 7.21*x[20] + 7.21*x[21] + 5.45*x[22] + 5.45*x[23] + 5.45*x[24] - 0.3*x[25] - 0.3*x[26] - 0.3*x[27] + 7.71*x[28] + 7.71*x[29] + 7.71*x[30] + 6.28*x[31] + 6.28*x[32] + 6.28*x[33] + 7.74*x[34] + 7.74*x[35] + 7.74*x[36] - 0.84*x[67] - 0.84*x[68] - 0.84*x[69] - 0.05*x[70] - 0.05*x[71] - 0.05*x[72] - 0.94*x[73] - 0.94*x[74] - 0.94*x[75] - 0.81*x[76] - 0.81*x[77] - 0.81*x[78] - 0.79*x[79] - 0.79*x[80] - 0.79*x[81] - 0.05*x[82] - 0.05*x[83] - 0.05*x[84] - 0.65*x[85] - 0.65*x[86] - 0.65*x[87] - 0.97*x[88] - 0.97*x[89] - 0.97*x[90] - 0.57*x[91] - 0.57*x[92] - 0.57*x[93] - 0.26*x[94] - 0.26*x[95] - 0.26*x[96] - 0.45*x[97] - 0.45*x[98] - 0.45*x[99] - 0.1*x[100] - 0.1*x[101] - 0.1*x[102])

    @NLconstraint(m, x[37]*x[55]-0.6*x[1]-0.2*x[13]+0.2*x[25]+0.2*x[28]+0.2*x[31] ==0.04)  #= e8: =#
    @NLconstraint(m, x[40]*x[58]-0.6*x[4]-0.2*x[16]-0.2*x[25]+0.7*x[34] ==0.07)  #= e9: =#
    @NLconstraint(m, x[43]*x[55]-0.4*x[1]-0.4*x[13]+0.5*x[25]+0.5*x[28]+0.5*x[31] ==0.1)  #= e10: =#
    @NLconstraint(m, x[46]*x[58]-0.4*x[4]-0.4*x[16]-0.5*x[25]+0.6*x[34] ==0.06)  #= e11: =#
    @NLconstraint(m, x[38]*x[56]-(x[37]*x[55]-(x[37]*x[26]+x[37]*x[29]+x[37]*x[32]))-0.6*x[2]-0.2*x[14] ==0)  #= e24: =#
    @NLconstraint(m, x[39]*x[57]-(x[38]*x[56]-(x[38]*x[27]+x[38]*x[30]+x[38]*x[33]))-0.6*x[3]-0.2*x[15] ==0)  #= e25: =#
    @NLconstraint(m, x[41]*x[59]-(x[40]*x[58]+x[37]*x[26]-x[40]*x[35])-0.6*x[5]-0.2*x[17] ==0)  #= e26: =#
    @NLconstraint(m, x[42]*x[60]-(x[41]*x[59]+x[38]*x[27]-x[41]*x[36])-0.6*x[6]-0.2*x[18] ==0)  #= e27: =#
    @NLconstraint(m, x[44]*x[56]-(x[43]*x[55]-(x[43]*x[26]+x[43]*x[29]+x[43]*x[32]))-0.4*x[2]-0.4*x[14] ==0)  #= e28: =#
    @NLconstraint(m, x[45]*x[57]-(x[44]*x[56]-(x[44]*x[27]+x[44]*x[30]+x[44]*x[33]))-0.4*x[3]-0.4*x[15] ==0)  #= e29: =#
    @NLconstraint(m, x[47]*x[59]-(x[46]*x[58]+x[43]*x[26]-x[46]*x[35])-0.4*x[5]-0.4*x[17] ==0)  #= e30: =#
    @NLconstraint(m, x[48]*x[60]-(x[47]*x[59]+x[44]*x[27]-x[47]*x[36])-0.4*x[6]-0.4*x[18] ==0)  #= e31: =#

    @constraint(m, x[1]+x[4]+x[7]+x[10]+x[49] ==1)  #= e2: =#
    @constraint(m, x[13]+x[16]+x[19]+x[22]+x[52] ==1.1)  #= e3: =#
    @constraint(m, -x[1]-x[13]+x[25]+x[28]+x[31]+x[55] ==0.2)  #= e4: =#
    @constraint(m, -x[4]-x[16]-x[25]+x[34]+x[58] ==0.1)  #= e5: =#
    @constraint(m, -x[7]-x[19]-x[28]-x[34]+x[61] ==1.55)  #= e6: =#
    @constraint(m, -x[10]-x[22]-x[31]+x[64] ==0.49)  #= e7: =#
    @constraint(m, x[2]+x[5]+x[8]+x[11]-x[49]+x[50] ==1)  #= e12: =#
    @constraint(m, x[3]+x[6]+x[9]+x[12]-x[50]+x[51] ==0)  #= e13: =#
    @constraint(m, x[14]+x[17]+x[20]+x[23]-x[52]+x[53] ==0.1)  #= e14: =#
    @constraint(m, x[15]+x[18]+x[21]+x[24]-x[53]+x[54] ==0.9)  #= e15: =#
    @constraint(m, -x[2]-x[14]+x[26]+x[29]+x[32]-x[55]+x[56] ==0)  #= e16: =#
    @constraint(m, -x[3]-x[15]+x[27]+x[30]+x[33]-x[56]+x[57] ==0)  #= e17: =#
    @constraint(m, -x[5]-x[17]-x[26]+x[35]-x[58]+x[59] ==0)  #= e18: =#
    @constraint(m, -x[6]-x[18]-x[27]+x[36]-x[59]+x[60] ==0)  #= e19: =#
    @constraint(m, -x[8]-x[20]-x[29]-x[35]-x[61]+x[62] ==-0.81)  #= e20: =#
    @constraint(m, -x[9]-x[21]-x[30]-x[36]-x[62]+x[63] ==-0.88)  #= e21: =#
    @constraint(m, -x[11]-x[23]-x[32]-x[64]+x[65] ==-0.14)  #= e22: =#
    @constraint(m, -x[12]-x[24]-x[33]-x[65]+x[66] ==-0.1)  #= e23: =#
    @constraint(m, x[1]-x[67]<=0)  #= e32: =#
    @constraint(m, x[2]-x[68]<=0)  #= e33: =#
    @constraint(m, x[3]-x[69]<=0)  #= e34: =#
    @constraint(m, x[4]-x[70]<=0)  #= e35: =#
    @constraint(m, x[5]-x[71]<=0)  #= e36: =#
    @constraint(m, x[6]-x[72]<=0)  #= e37: =#
    @constraint(m, x[7]-x[73]<=0)  #= e38: =#
    @constraint(m, x[8]-x[74]<=0)  #= e39: =#
    @constraint(m, x[9]-x[75]<=0)  #= e40: =#
    @constraint(m, x[10]-x[76]<=0)  #= e41: =#
    @constraint(m, x[11]-x[77]<=0)  #= e42: =#
    @constraint(m, x[12]-x[78]<=0)  #= e43: =#
    @constraint(m, x[13]-x[79]<=0)  #= e44: =#
    @constraint(m, x[14]-x[80]<=0)  #= e45: =#
    @constraint(m, x[15]-x[81]<=0)  #= e46: =#
    @constraint(m, x[16]-x[82]<=0)  #= e47: =#
    @constraint(m, x[17]-x[83]<=0)  #= e48: =#
    @constraint(m, x[18]-x[84]<=0)  #= e49: =#
    @constraint(m, x[19]-x[85]<=0)  #= e50: =#
    @constraint(m, x[20]-x[86]<=0)  #= e51: =#
    @constraint(m, x[21]-x[87]<=0)  #= e52: =#
    @constraint(m, x[22]-x[88]<=0)  #= e53: =#
    @constraint(m, x[23]-x[89]<=0)  #= e54: =#
    @constraint(m, x[24]-x[90]<=0)  #= e55: =#
    @constraint(m, x[25]-x[91]<=0)  #= e56: =#
    @constraint(m, x[26]-x[92]<=0)  #= e57: =#
    @constraint(m, x[27]-x[93]<=0)  #= e58: =#
    @constraint(m, x[28]-x[94]<=0)  #= e59: =#
    @constraint(m, x[29]-x[95]<=0)  #= e60: =#
    @constraint(m, x[30]-x[96]<=0)  #= e61: =#
    @constraint(m, x[31]-x[97]<=0)  #= e62: =#
    @constraint(m, x[32]-x[98]<=0)  #= e63: =#
    @constraint(m, x[33]-x[99]<=0)  #= e64: =#
    @constraint(m, x[34]-x[100]<=0)  #= e65: =#
    @constraint(m, x[35]-x[101]<=0)  #= e66: =#
    @constraint(m, x[36]-x[102]<=0)  #= e67: =#
    @constraint(m, x[1]>=0)  #= e68: =#
    @constraint(m, x[2]>=0)  #= e69: =#
    @constraint(m, x[3]>=0)  #= e70: =#
    @constraint(m, x[4]>=0)  #= e71: =#
    @constraint(m, x[5]>=0)  #= e72: =#
    @constraint(m, x[6]>=0)  #= e73: =#
    @constraint(m, x[7]>=0)  #= e74: =#
    @constraint(m, x[8]>=0)  #= e75: =#
    @constraint(m, x[9]>=0)  #= e76: =#
    @constraint(m, x[10]>=0)  #= e77: =#
    @constraint(m, x[11]>=0)  #= e78: =#
    @constraint(m, x[12]>=0)  #= e79: =#
    @constraint(m, x[13]>=0)  #= e80: =#
    @constraint(m, x[14]>=0)  #= e81: =#
    @constraint(m, x[15]>=0)  #= e82: =#
    @constraint(m, x[16]>=0)  #= e83: =#
    @constraint(m, x[17]>=0)  #= e84: =#
    @constraint(m, x[18]>=0)  #= e85: =#
    @constraint(m, x[19]>=0)  #= e86: =#
    @constraint(m, x[20]>=0)  #= e87: =#
    @constraint(m, x[21]>=0)  #= e88: =#
    @constraint(m, x[22]>=0)  #= e89: =#
    @constraint(m, x[23]>=0)  #= e90: =#
    @constraint(m, x[24]>=0)  #= e91: =#
    @constraint(m, x[25]>=0)  #= e92: =#
    @constraint(m, x[26]>=0)  #= e93: =#
    @constraint(m, x[27]>=0)  #= e94: =#
    @constraint(m, x[28]>=0)  #= e95: =#
    @constraint(m, x[29]>=0)  #= e96: =#
    @constraint(m, x[30]>=0)  #= e97: =#
    @constraint(m, x[31]>=0)  #= e98: =#
    @constraint(m, x[32]>=0)  #= e99: =#
    @constraint(m, x[33]>=0)  #= e100: =#
    @constraint(m, x[34]>=0)  #= e101: =#
    @constraint(m, x[35]>=0)  #= e102: =#
    @constraint(m, x[36]>=0)  #= e103: =#
    @constraint(m, x[73]<=1.5)  #= e104: =#
    @constraint(m, x[74]<=1.5)  #= e105: =#
    @constraint(m, x[75]<=1.5)  #= e106: =#
    @constraint(m, x[76]<=0.6)  #= e107: =#
    @constraint(m, x[77]<=0.6)  #= e108: =#
    @constraint(m, x[78]<=0.6)  #= e109: =#
    @constraint(m, x[85]<=1.1)  #= e110: =#
    @constraint(m, x[86]<=1.1)  #= e111: =#
    @constraint(m, x[87]<=1.1)  #= e112: =#
    @constraint(m, x[88]<=0.2)  #= e113: =#
    @constraint(m, x[89]<=0.2)  #= e114: =#
    @constraint(m, x[90]<=0.2)  #= e115: =#
    @constraint(m, x[73]<=1)  #= e116: =#
    @constraint(m, x[74]<=1)  #= e117: =#
    @constraint(m, x[75]<=1)  #= e118: =#
    @constraint(m, x[76]<=0.8)  #= e119: =#
    @constraint(m, x[77]<=0.8)  #= e120: =#
    @constraint(m, x[78]<=0.8)  #= e121: =#
    @constraint(m, x[85]<=1)  #= e122: =#
    @constraint(m, x[86]<=1)  #= e123: =#
    @constraint(m, x[87]<=1)  #= e124: =#
    @constraint(m, x[88]<=0.8)  #= e125: =#
    @constraint(m, x[89]<=0.8)  #= e126: =#
    @constraint(m, x[90]<=0.8)  #= e127: =#
    @constraint(m, -x[73]>=-1.3)  #= e128: =#
    @constraint(m, -x[74]>=-1.3)  #= e129: =#
    @constraint(m, -x[75]>=-1.3)  #= e130: =#
    @constraint(m, -x[76]>=-1.4)  #= e131: =#
    @constraint(m, -x[77]>=-1.4)  #= e132: =#
    @constraint(m, -x[78]>=-1.4)  #= e133: =#
    @constraint(m, -x[85]>=-1.7)  #= e134: =#
    @constraint(m, -x[86]>=-1.7)  #= e135: =#
    @constraint(m, -x[87]>=-1.7)  #= e136: =#
    @constraint(m, -x[88]>=-1.8)  #= e137: =#
    @constraint(m, -x[89]>=-1.8)  #= e138: =#
    @constraint(m, -x[90]>=-1.8)  #= e139: =#
    @constraint(m, -x[73]>=-1)  #= e140: =#
    @constraint(m, -x[74]>=-1)  #= e141: =#
    @constraint(m, -x[75]>=-1)  #= e142: =#
    @constraint(m, -x[76]>=-1.4)  #= e143: =#
    @constraint(m, -x[77]>=-1.4)  #= e144: =#
    @constraint(m, -x[78]>=-1.4)  #= e145: =#
    @constraint(m, -x[85]>=-1)  #= e146: =#
    @constraint(m, -x[86]>=-1)  #= e147: =#
    @constraint(m, -x[87]>=-1)  #= e148: =#
    @constraint(m, -x[88]>=-1.4)  #= e149: =#
    @constraint(m, -x[89]>=-1.4)  #= e150: =#
    @constraint(m, -x[90]>=-1.4)  #= e151: =#
    @constraint(m, -x[37]+x[95]<=0.9)  #= e152: =#
    @constraint(m, -x[38]+x[96]<=0.9)  #= e153: =#
    @constraint(m, -x[37]+x[98]<=0)  #= e154: =#
    @constraint(m, -x[38]+x[99]<=0)  #= e155: =#
    @constraint(m, -x[40]+x[101]<=0.9)  #= e156: =#
    @constraint(m, -x[41]+x[102]<=0.9)  #= e157: =#
    @constraint(m, -x[43]+x[95]<=0.6)  #= e158: =#
    @constraint(m, -x[44]+x[96]<=0.6)  #= e159: =#
    @constraint(m, -x[43]+x[98]<=0.4)  #= e160: =#
    @constraint(m, -x[44]+x[99]<=0.4)  #= e161: =#
    @constraint(m, -x[46]+x[101]<=0.6)  #= e162: =#
    @constraint(m, -x[47]+x[102]<=0.6)  #= e163: =#
    @constraint(m, -x[37]-x[95]>=-1.9)  #= e164: =#
    @constraint(m, -x[38]-x[96]>=-1.9)  #= e165: =#
    @constraint(m, -x[37]-x[98]>=-2)  #= e166: =#
    @constraint(m, -x[38]-x[99]>=-2)  #= e167: =#
    @constraint(m, -x[40]-x[101]>=-1.9)  #= e168: =#
    @constraint(m, -x[41]-x[102]>=-1.9)  #= e169: =#
    @constraint(m, -x[43]-x[95]>=-1.4)  #= e170: =#
    @constraint(m, -x[44]-x[96]>=-1.4)  #= e171: =#
    @constraint(m, -x[43]-x[98]>=-1.8)  #= e172: =#
    @constraint(m, -x[44]-x[99]>=-1.8)  #= e173: =#
    @constraint(m, -x[46]-x[101]>=-1.4)  #= e174: =#
    @constraint(m, -x[47]-x[102]>=-1.4)  #= e175: =#
    @constraint(m, x[94]<=1.1)  #= e176: =#
    @constraint(m, x[97]<=0.2)  #= e177: =#
    @constraint(m, x[100]<=1.6)  #= e178: =#
    @constraint(m, x[94]<=1.1)  #= e179: =#
    @constraint(m, x[97]<=0.9)  #= e180: =#
    @constraint(m, x[100]<=1.2)  #= e181: =#
    @constraint(m, -x[94]>=-1.7)  #= e182: =#
    @constraint(m, -x[97]>=-1.8)  #= e183: =#
    @constraint(m, -x[100]>=-1.2)  #= e184: =#
    @constraint(m, -x[94]>=-0.9)  #= e185: =#
    @constraint(m, -x[97]>=-1.3)  #= e186: =#
    @constraint(m, -x[100]>=-0.8)  #= e187: =#
    @constraint(m, x[67]+x[91]<=1)  #= e188: =#
    @constraint(m, x[68]+x[92]<=1)  #= e189: =#
    @constraint(m, x[69]+x[93]<=1)  #= e190: =#
    @constraint(m, x[67]+x[94]<=1)  #= e191: =#
    @constraint(m, x[68]+x[95]<=1)  #= e192: =#
    @constraint(m, x[69]+x[96]<=1)  #= e193: =#
    @constraint(m, x[67]+x[97]<=1)  #= e194: =#
    @constraint(m, x[68]+x[98]<=1)  #= e195: =#
    @constraint(m, x[69]+x[99]<=1)  #= e196: =#
    @constraint(m, x[79]+x[91]<=1)  #= e197: =#
    @constraint(m, x[80]+x[92]<=1)  #= e198: =#
    @constraint(m, x[81]+x[93]<=1)  #= e199: =#
    @constraint(m, x[79]+x[94]<=1)  #= e200: =#
    @constraint(m, x[80]+x[95]<=1)  #= e201: =#
    @constraint(m, x[81]+x[96]<=1)  #= e202: =#
    @constraint(m, x[79]+x[97]<=1)  #= e203: =#
    @constraint(m, x[80]+x[98]<=1)  #= e204: =#
    @constraint(m, x[81]+x[99]<=1)  #= e205: =#
    @constraint(m, x[70]+x[100]<=1)  #= e206: =#
    @constraint(m, x[71]+x[101]<=1)  #= e207: =#
    @constraint(m, x[72]+x[102]<=1)  #= e208: =#
    @constraint(m, x[82]+x[100]<=1)  #= e209: =#
    @constraint(m, x[83]+x[101]<=1)  #= e210: =#
    @constraint(m, x[84]+x[102]<=1)  #= e211: =#
    @constraint(m, x[91]+x[100]<=1)  #= e212: =#
    @constraint(m, x[92]+x[101]<=1)  #= e213: =#
    @constraint(m, x[93]+x[102]<=1)  #= e214: =#

    return m
end

function blend029_gl(;solver=nothing)

    m = Model(solver)

    @variable(m, x[1:102])
    for i=67:102
      JuMP.set_binary(x[i])
    end

    for i=1:48
      set_lower_bound(x[i], 0)
      set_upper_bound(x[i], 1)
    end
    for i=49:66
      set_lower_bound(x[i], 0)
      set_upper_bound(x[i], 2)
    end

    @objective(m, Max, - 1.74*x[1] - 1.74*x[2] - 1.74*x[3] - 1.45*x[4] - 1.45*x[5] - 1.45*x[6] + 7.38*x[7] + 7.38*x[8] + 7.38*x[9] + 5.6*x[10] + 5.6*x[11] + 5.6*x[12] - 1.7*x[13] - 1.7*x[14] - 1.7*x[15] - 1.18*x[16] - 1.18*x[17] - 1.18*x[18] + 7.21*x[19] + 7.21*x[20] + 7.21*x[21] + 5.45*x[22] + 5.45*x[23] + 5.45*x[24] - 0.3*x[25] - 0.3*x[26] - 0.3*x[27] + 7.71*x[28] + 7.71*x[29] + 7.71*x[30] + 6.28*x[31] + 6.28*x[32] + 6.28*x[33] + 7.74*x[34] + 7.74*x[35] + 7.74*x[36] - 0.84*x[67] - 0.84*x[68] - 0.84*x[69] - 0.05*x[70] - 0.05*x[71] - 0.05*x[72] - 0.94*x[73] - 0.94*x[74] - 0.94*x[75] - 0.81*x[76] - 0.81*x[77] - 0.81*x[78] - 0.79*x[79] - 0.79*x[80] - 0.79*x[81] - 0.05*x[82] - 0.05*x[83] - 0.05*x[84] - 0.65*x[85] - 0.65*x[86] - 0.65*x[87] - 0.97*x[88] - 0.97*x[89] - 0.97*x[90] - 0.57*x[91] - 0.57*x[92] - 0.57*x[93] - 0.26*x[94] - 0.26*x[95] - 0.26*x[96] - 0.45*x[97] - 0.45*x[98] - 0.45*x[99] - 0.1*x[100] - 0.1*x[101] - 0.1*x[102])

    @NLconstraint(m, x[37]*x[55]-0.6*x[1]-0.2*x[13]+0.2*x[25]+0.2*x[28]+0.2*x[31] >=0.04)  #= e8: =#
    @NLconstraint(m, x[40]*x[58]-0.6*x[4]-0.2*x[16]-0.2*x[25]+0.7*x[34] >=0.07)  #= e9: =#
    @NLconstraint(m, x[43]*x[55]-0.4*x[1]-0.4*x[13]+0.5*x[25]+0.5*x[28]+0.5*x[31] >=0.1)  #= e10: =#
    @NLconstraint(m, x[46]*x[58]-0.4*x[4]-0.4*x[16]-0.5*x[25]+0.6*x[34] >=0.06)  #= e11: =#
    @NLconstraint(m, x[38]*x[56]-(x[37]*x[55]-(x[37]*x[26]+x[37]*x[29]+x[37]*x[32]))-0.6*x[2]-0.2*x[14] >=0)  #= e24: =#
    @NLconstraint(m, x[39]*x[57]-(x[38]*x[56]-(x[38]*x[27]+x[38]*x[30]+x[38]*x[33]))-0.6*x[3]-0.2*x[15] >=0)  #= e25: =#
    @NLconstraint(m, x[41]*x[59]-(x[40]*x[58]+x[37]*x[26]-x[40]*x[35])-0.6*x[5]-0.2*x[17] >=0)  #= e26: =#
    @NLconstraint(m, x[42]*x[60]-(x[41]*x[59]+x[38]*x[27]-x[41]*x[36])-0.6*x[6]-0.2*x[18] >=0)  #= e27: =#
    @NLconstraint(m, x[44]*x[56]-(x[43]*x[55]-(x[43]*x[26]+x[43]*x[29]+x[43]*x[32]))-0.4*x[2]-0.4*x[14] >=0)  #= e28: =#
    @NLconstraint(m, x[45]*x[57]-(x[44]*x[56]-(x[44]*x[27]+x[44]*x[30]+x[44]*x[33]))-0.4*x[3]-0.4*x[15] >=0)  #= e29: =#
    @NLconstraint(m, x[47]*x[59]-(x[46]*x[58]+x[43]*x[26]-x[46]*x[35])-0.4*x[5]-0.4*x[17] >=0)  #= e30: =#
    @NLconstraint(m, x[48]*x[60]-(x[47]*x[59]+x[44]*x[27]-x[47]*x[36])-0.4*x[6]-0.4*x[18] >=0)  #= e31: =#

    @NLconstraint(m, x[37]*x[55]-0.6*x[1]-0.2*x[13]+0.2*x[25]+0.2*x[28]+0.2*x[31] <=0.04)  #= e8: =#
    @NLconstraint(m, x[40]*x[58]-0.6*x[4]-0.2*x[16]-0.2*x[25]+0.7*x[34] <=0.07)  #= e9: =#
    @NLconstraint(m, x[43]*x[55]-0.4*x[1]-0.4*x[13]+0.5*x[25]+0.5*x[28]+0.5*x[31] <=0.1)  #= e10: =#
    @NLconstraint(m, x[46]*x[58]-0.4*x[4]-0.4*x[16]-0.5*x[25]+0.6*x[34] <=0.06)  #= e11: =#
    @NLconstraint(m, x[38]*x[56]-(x[37]*x[55]-(x[37]*x[26]+x[37]*x[29]+x[37]*x[32]))-0.6*x[2]-0.2*x[14] <=0)  #= e24: =#
    @NLconstraint(m, x[39]*x[57]-(x[38]*x[56]-(x[38]*x[27]+x[38]*x[30]+x[38]*x[33]))-0.6*x[3]-0.2*x[15] <=0)  #= e25: =#
    @NLconstraint(m, x[41]*x[59]-(x[40]*x[58]+x[37]*x[26]-x[40]*x[35])-0.6*x[5]-0.2*x[17] <=0)  #= e26: =#
    @NLconstraint(m, x[42]*x[60]-(x[41]*x[59]+x[38]*x[27]-x[41]*x[36])-0.6*x[6]-0.2*x[18] <=0)  #= e27: =#
    @NLconstraint(m, x[44]*x[56]-(x[43]*x[55]-(x[43]*x[26]+x[43]*x[29]+x[43]*x[32]))-0.4*x[2]-0.4*x[14] <=0)  #= e28: =#
    @NLconstraint(m, x[45]*x[57]-(x[44]*x[56]-(x[44]*x[27]+x[44]*x[30]+x[44]*x[33]))-0.4*x[3]-0.4*x[15] <=0)  #= e29: =#
    @NLconstraint(m, x[47]*x[59]-(x[46]*x[58]+x[43]*x[26]-x[46]*x[35])-0.4*x[5]-0.4*x[17] <=0)  #= e30: =#
    @NLconstraint(m, x[48]*x[60]-(x[47]*x[59]+x[44]*x[27]-x[47]*x[36])-0.4*x[6]-0.4*x[18] <=0)  #= e31: =#

    warmstarter = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.53351, 0.60649, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 1.0, 2.0, 1.0, 1.1, 0.66649, 0.96, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 1.55, 1.27351, 2.0, 0.49, 0.35, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for i in 1:102
      set_start_value(x[i], warmstarter[i])
    end

    @constraint(m, x[1]+x[4]+x[7]+x[10]+x[49] ==1)  #= e2: =#
    @constraint(m, x[13]+x[16]+x[19]+x[22]+x[52] ==1.1)  #= e3: =#
    @constraint(m, -x[1]-x[13]+x[25]+x[28]+x[31]+x[55] ==0.2)  #= e4: =#
    @constraint(m, -x[4]-x[16]-x[25]+x[34]+x[58] ==0.1)  #= e5: =#
    @constraint(m, -x[7]-x[19]-x[28]-x[34]+x[61] ==1.55)  #= e6: =#
    @constraint(m, -x[10]-x[22]-x[31]+x[64] ==0.49)  #= e7: =#
    @constraint(m, x[2]+x[5]+x[8]+x[11]-x[49]+x[50] ==1)  #= e12: =#
    @constraint(m, x[3]+x[6]+x[9]+x[12]-x[50]+x[51] ==0)  #= e13: =#
    @constraint(m, x[14]+x[17]+x[20]+x[23]-x[52]+x[53] ==0.1)  #= e14: =#
    @constraint(m, x[15]+x[18]+x[21]+x[24]-x[53]+x[54] ==0.9)  #= e15: =#
    @constraint(m, -x[2]-x[14]+x[26]+x[29]+x[32]-x[55]+x[56] ==0)  #= e16: =#
    @constraint(m, -x[3]-x[15]+x[27]+x[30]+x[33]-x[56]+x[57] ==0)  #= e17: =#
    @constraint(m, -x[5]-x[17]-x[26]+x[35]-x[58]+x[59] ==0)  #= e18: =#
    @constraint(m, -x[6]-x[18]-x[27]+x[36]-x[59]+x[60] ==0)  #= e19: =#
    @constraint(m, -x[8]-x[20]-x[29]-x[35]-x[61]+x[62] ==-0.81)  #= e20: =#
    @constraint(m, -x[9]-x[21]-x[30]-x[36]-x[62]+x[63] ==-0.88)  #= e21: =#
    @constraint(m, -x[11]-x[23]-x[32]-x[64]+x[65] ==-0.14)  #= e22: =#
    @constraint(m, -x[12]-x[24]-x[33]-x[65]+x[66] ==-0.1)  #= e23: =#
    @constraint(m, x[1]-x[67]<=0)  #= e32: =#
    @constraint(m, x[2]-x[68]<=0)  #= e33: =#
    @constraint(m, x[3]-x[69]<=0)  #= e34: =#
    @constraint(m, x[4]-x[70]<=0)  #= e35: =#
    @constraint(m, x[5]-x[71]<=0)  #= e36: =#
    @constraint(m, x[6]-x[72]<=0)  #= e37: =#
    @constraint(m, x[7]-x[73]<=0)  #= e38: =#
    @constraint(m, x[8]-x[74]<=0)  #= e39: =#
    @constraint(m, x[9]-x[75]<=0)  #= e40: =#
    @constraint(m, x[10]-x[76]<=0)  #= e41: =#
    @constraint(m, x[11]-x[77]<=0)  #= e42: =#
    @constraint(m, x[12]-x[78]<=0)  #= e43: =#
    @constraint(m, x[13]-x[79]<=0)  #= e44: =#
    @constraint(m, x[14]-x[80]<=0)  #= e45: =#
    @constraint(m, x[15]-x[81]<=0)  #= e46: =#
    @constraint(m, x[16]-x[82]<=0)  #= e47: =#
    @constraint(m, x[17]-x[83]<=0)  #= e48: =#
    @constraint(m, x[18]-x[84]<=0)  #= e49: =#
    @constraint(m, x[19]-x[85]<=0)  #= e50: =#
    @constraint(m, x[20]-x[86]<=0)  #= e51: =#
    @constraint(m, x[21]-x[87]<=0)  #= e52: =#
    @constraint(m, x[22]-x[88]<=0)  #= e53: =#
    @constraint(m, x[23]-x[89]<=0)  #= e54: =#
    @constraint(m, x[24]-x[90]<=0)  #= e55: =#
    @constraint(m, x[25]-x[91]<=0)  #= e56: =#
    @constraint(m, x[26]-x[92]<=0)  #= e57: =#
    @constraint(m, x[27]-x[93]<=0)  #= e58: =#
    @constraint(m, x[28]-x[94]<=0)  #= e59: =#
    @constraint(m, x[29]-x[95]<=0)  #= e60: =#
    @constraint(m, x[30]-x[96]<=0)  #= e61: =#
    @constraint(m, x[31]-x[97]<=0)  #= e62: =#
    @constraint(m, x[32]-x[98]<=0)  #= e63: =#
    @constraint(m, x[33]-x[99]<=0)  #= e64: =#
    @constraint(m, x[34]-x[100]<=0)  #= e65: =#
    @constraint(m, x[35]-x[101]<=0)  #= e66: =#
    @constraint(m, x[36]-x[102]<=0)  #= e67: =#
    @constraint(m, x[1]>=0)  #= e68: =#
    @constraint(m, x[2]>=0)  #= e69: =#
    @constraint(m, x[3]>=0)  #= e70: =#
    @constraint(m, x[4]>=0)  #= e71: =#
    @constraint(m, x[5]>=0)  #= e72: =#
    @constraint(m, x[6]>=0)  #= e73: =#
    @constraint(m, x[7]>=0)  #= e74: =#
    @constraint(m, x[8]>=0)  #= e75: =#
    @constraint(m, x[9]>=0)  #= e76: =#
    @constraint(m, x[10]>=0)  #= e77: =#
    @constraint(m, x[11]>=0)  #= e78: =#
    @constraint(m, x[12]>=0)  #= e79: =#
    @constraint(m, x[13]>=0)  #= e80: =#
    @constraint(m, x[14]>=0)  #= e81: =#
    @constraint(m, x[15]>=0)  #= e82: =#
    @constraint(m, x[16]>=0)  #= e83: =#
    @constraint(m, x[17]>=0)  #= e84: =#
    @constraint(m, x[18]>=0)  #= e85: =#
    @constraint(m, x[19]>=0)  #= e86: =#
    @constraint(m, x[20]>=0)  #= e87: =#
    @constraint(m, x[21]>=0)  #= e88: =#
    @constraint(m, x[22]>=0)  #= e89: =#
    @constraint(m, x[23]>=0)  #= e90: =#
    @constraint(m, x[24]>=0)  #= e91: =#
    @constraint(m, x[25]>=0)  #= e92: =#
    @constraint(m, x[26]>=0)  #= e93: =#
    @constraint(m, x[27]>=0)  #= e94: =#
    @constraint(m, x[28]>=0)  #= e95: =#
    @constraint(m, x[29]>=0)  #= e96: =#
    @constraint(m, x[30]>=0)  #= e97: =#
    @constraint(m, x[31]>=0)  #= e98: =#
    @constraint(m, x[32]>=0)  #= e99: =#
    @constraint(m, x[33]>=0)  #= e100: =#
    @constraint(m, x[34]>=0)  #= e101: =#
    @constraint(m, x[35]>=0)  #= e102: =#
    @constraint(m, x[36]>=0)  #= e103: =#
    @constraint(m, x[73]<=1.5)  #= e104: =#
    @constraint(m, x[74]<=1.5)  #= e105: =#
    @constraint(m, x[75]<=1.5)  #= e106: =#
    @constraint(m, x[76]<=0.6)  #= e107: =#
    @constraint(m, x[77]<=0.6)  #= e108: =#
    @constraint(m, x[78]<=0.6)  #= e109: =#
    @constraint(m, x[85]<=1.1)  #= e110: =#
    @constraint(m, x[86]<=1.1)  #= e111: =#
    @constraint(m, x[87]<=1.1)  #= e112: =#
    @constraint(m, x[88]<=0.2)  #= e113: =#
    @constraint(m, x[89]<=0.2)  #= e114: =#
    @constraint(m, x[90]<=0.2)  #= e115: =#
    @constraint(m, x[73]<=1)  #= e116: =#
    @constraint(m, x[74]<=1)  #= e117: =#
    @constraint(m, x[75]<=1)  #= e118: =#
    @constraint(m, x[76]<=0.8)  #= e119: =#
    @constraint(m, x[77]<=0.8)  #= e120: =#
    @constraint(m, x[78]<=0.8)  #= e121: =#
    @constraint(m, x[85]<=1)  #= e122: =#
    @constraint(m, x[86]<=1)  #= e123: =#
    @constraint(m, x[87]<=1)  #= e124: =#
    @constraint(m, x[88]<=0.8)  #= e125: =#
    @constraint(m, x[89]<=0.8)  #= e126: =#
    @constraint(m, x[90]<=0.8)  #= e127: =#
    @constraint(m, -x[73]>=-1.3)  #= e128: =#
    @constraint(m, -x[74]>=-1.3)  #= e129: =#
    @constraint(m, -x[75]>=-1.3)  #= e130: =#
    @constraint(m, -x[76]>=-1.4)  #= e131: =#
    @constraint(m, -x[77]>=-1.4)  #= e132: =#
    @constraint(m, -x[78]>=-1.4)  #= e133: =#
    @constraint(m, -x[85]>=-1.7)  #= e134: =#
    @constraint(m, -x[86]>=-1.7)  #= e135: =#
    @constraint(m, -x[87]>=-1.7)  #= e136: =#
    @constraint(m, -x[88]>=-1.8)  #= e137: =#
    @constraint(m, -x[89]>=-1.8)  #= e138: =#
    @constraint(m, -x[90]>=-1.8)  #= e139: =#
    @constraint(m, -x[73]>=-1)  #= e140: =#
    @constraint(m, -x[74]>=-1)  #= e141: =#
    @constraint(m, -x[75]>=-1)  #= e142: =#
    @constraint(m, -x[76]>=-1.4)  #= e143: =#
    @constraint(m, -x[77]>=-1.4)  #= e144: =#
    @constraint(m, -x[78]>=-1.4)  #= e145: =#
    @constraint(m, -x[85]>=-1)  #= e146: =#
    @constraint(m, -x[86]>=-1)  #= e147: =#
    @constraint(m, -x[87]>=-1)  #= e148: =#
    @constraint(m, -x[88]>=-1.4)  #= e149: =#
    @constraint(m, -x[89]>=-1.4)  #= e150: =#
    @constraint(m, -x[90]>=-1.4)  #= e151: =#
    @constraint(m, -x[37]+x[95]<=0.9)  #= e152: =#
    @constraint(m, -x[38]+x[96]<=0.9)  #= e153: =#
    @constraint(m, -x[37]+x[98]<=0)  #= e154: =#
    @constraint(m, -x[38]+x[99]<=0)  #= e155: =#
    @constraint(m, -x[40]+x[101]<=0.9)  #= e156: =#
    @constraint(m, -x[41]+x[102]<=0.9)  #= e157: =#
    @constraint(m, -x[43]+x[95]<=0.6)  #= e158: =#
    @constraint(m, -x[44]+x[96]<=0.6)  #= e159: =#
    @constraint(m, -x[43]+x[98]<=0.4)  #= e160: =#
    @constraint(m, -x[44]+x[99]<=0.4)  #= e161: =#
    @constraint(m, -x[46]+x[101]<=0.6)  #= e162: =#
    @constraint(m, -x[47]+x[102]<=0.6)  #= e163: =#
    @constraint(m, -x[37]-x[95]>=-1.9)  #= e164: =#
    @constraint(m, -x[38]-x[96]>=-1.9)  #= e165: =#
    @constraint(m, -x[37]-x[98]>=-2)  #= e166: =#
    @constraint(m, -x[38]-x[99]>=-2)  #= e167: =#
    @constraint(m, -x[40]-x[101]>=-1.9)  #= e168: =#
    @constraint(m, -x[41]-x[102]>=-1.9)  #= e169: =#
    @constraint(m, -x[43]-x[95]>=-1.4)  #= e170: =#
    @constraint(m, -x[44]-x[96]>=-1.4)  #= e171: =#
    @constraint(m, -x[43]-x[98]>=-1.8)  #= e172: =#
    @constraint(m, -x[44]-x[99]>=-1.8)  #= e173: =#
    @constraint(m, -x[46]-x[101]>=-1.4)  #= e174: =#
    @constraint(m, -x[47]-x[102]>=-1.4)  #= e175: =#
    @constraint(m, x[94]<=1.1)  #= e176: =#
    @constraint(m, x[97]<=0.2)  #= e177: =#
    @constraint(m, x[100]<=1.6)  #= e178: =#
    @constraint(m, x[94]<=1.1)  #= e179: =#
    @constraint(m, x[97]<=0.9)  #= e180: =#
    @constraint(m, x[100]<=1.2)  #= e181: =#
    @constraint(m, -x[94]>=-1.7)  #= e182: =#
    @constraint(m, -x[97]>=-1.8)  #= e183: =#
    @constraint(m, -x[100]>=-1.2)  #= e184: =#
    @constraint(m, -x[94]>=-0.9)  #= e185: =#
    @constraint(m, -x[97]>=-1.3)  #= e186: =#
    @constraint(m, -x[100]>=-0.8)  #= e187: =#
    @constraint(m, x[67]+x[91]<=1)  #= e188: =#
    @constraint(m, x[68]+x[92]<=1)  #= e189: =#
    @constraint(m, x[69]+x[93]<=1)  #= e190: =#
    @constraint(m, x[67]+x[94]<=1)  #= e191: =#
    @constraint(m, x[68]+x[95]<=1)  #= e192: =#
    @constraint(m, x[69]+x[96]<=1)  #= e193: =#
    @constraint(m, x[67]+x[97]<=1)  #= e194: =#
    @constraint(m, x[68]+x[98]<=1)  #= e195: =#
    @constraint(m, x[69]+x[99]<=1)  #= e196: =#
    @constraint(m, x[79]+x[91]<=1)  #= e197: =#
    @constraint(m, x[80]+x[92]<=1)  #= e198: =#
    @constraint(m, x[81]+x[93]<=1)  #= e199: =#
    @constraint(m, x[79]+x[94]<=1)  #= e200: =#
    @constraint(m, x[80]+x[95]<=1)  #= e201: =#
    @constraint(m, x[81]+x[96]<=1)  #= e202: =#
    @constraint(m, x[79]+x[97]<=1)  #= e203: =#
    @constraint(m, x[80]+x[98]<=1)  #= e204: =#
    @constraint(m, x[81]+x[99]<=1)  #= e205: =#
    @constraint(m, x[70]+x[100]<=1)  #= e206: =#
    @constraint(m, x[71]+x[101]<=1)  #= e207: =#
    @constraint(m, x[72]+x[102]<=1)  #= e208: =#
    @constraint(m, x[82]+x[100]<=1)  #= e209: =#
    @constraint(m, x[83]+x[101]<=1)  #= e210: =#
    @constraint(m, x[84]+x[102]<=1)  #= e211: =#
    @constraint(m, x[91]+x[100]<=1)  #= e212: =#
    @constraint(m, x[92]+x[101]<=1)  #= e213: =#
    @constraint(m, x[93]+x[102]<=1)  #= e214: =#

    return m
end
