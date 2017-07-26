using POD, JuMP, Gurobi, Ipopt, MathProgBase, AmplNLWriter, CoinOptServices

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
