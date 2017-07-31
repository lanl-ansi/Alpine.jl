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
