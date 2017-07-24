using POD, JuMP, Gurobi, AmplNLWriter, CoinOptServices, MathProgBase

function meanvarx(;verbose=false, solver=nothing)

    if solver==nothing
        m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(["bonmin.nlp_log_level=0"; "bonmin.bb_log_level=0"]),
                                   mip_solver=GurobiSolver(OutputFlag=1),
                                   rel_gap=0.001, log_level=100,
                                   discretization_ratio=4,
                                   bilinear_convexhull=false,
                                   presolve_bt_output_tol=1e-2,
                                   presolve_bound_tightening=false))
    else
        m = Model(solver=solver)
    end

    @variable(m, x[1:35]>=0)
    for i=22:35
        setcategory(x[i], :Bin)
    end

    # Artifical upper bounds
    for i=1:7
        setupperbound(x[i], 1.0)
    end
    # Artifical upper bounds

    @NLobjective(m, Min, 42.18*x[1]*x[1] + 40.36*x[1]*x[2] + 21.76*x[1]*x[3] + 10.6*x[1]*x[4] + 24.64*x[1]*x[5] + 47.68*x[1]*x[6] + 34.82*x[1]*x[7] + 70.89*x[2]*x[2] + 43.16*x[2]*x[3] + 30.82*x[2]*x[4] + 46.48*x[2]*x[5] + 47.6*x[2]*x[6] + 25.24*x[2]*x[7] + 25.51*x[3]*x[3] + 19.2*x[3]*x[4] + 45.26*x[3]*x[5] + 26.44*x[3]*x[6] + 9.4*x[3]*x[7] + 22.33*x[4]*x[4] + 20.64*x[4]*x[5] + 20.92*x[4]*x[6] + 2*x[4]*x[7] + 30.01*x[5]*x[5] + 32.72*x[5]*x[6] + 14.4*x[5]*x[7] + 42.23*x[6]*x[6] + 19.8*x[6]*x[7] + 16.42*x[7]*x[7] - 0.06435*x[1] - 0.0548*x[2] - 0.02505*x[3] - 0.0762*x[4] - 0.03815*x[5] - 0.0927*x[6] - 0.031*x[7])

    @constraint(m, x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]==1)  #= e1: =#
    @constraint(m, x[1]-x[8]+x[15] ==0.2)  #= e2: =#
    @constraint(m, x[2]-x[9]+x[16] ==0.2)  #= e3: =#
    @constraint(m, x[3]-x[10]+x[17] ==0)  #= e4: =#
    @constraint(m, x[4]-x[11]+x[18] ==0)  #= e5: =#
    @constraint(m, x[5]-x[12]+x[19] ==0.2)  #= e6: =#
    @constraint(m, x[6]-x[13]+x[20] ==0.2)  #= e7: =#
    @constraint(m, x[7]-x[14]+x[21] ==0.2)  #= e8: =#
    @constraint(m, x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]<=0.3)  #= e9: =#
    @constraint(m, x[8]-0.11*x[22]<=0)  #= e10: =#
    @constraint(m, x[9]-0.1*x[23]<=0)  #= e11: =#
    @constraint(m, x[10]-0.07*x[24]<=0)  #= e12: =#
    @constraint(m, x[11]-0.11*x[25]<=0)  #= e13: =#
    @constraint(m, x[12]-0.2*x[26]<=0)  #= e14: =#
    @constraint(m, x[13]-0.1*x[27]<=0)  #= e15: =#
    @constraint(m, x[14]-0.1*x[28]<=0)  #= e16: =#
    @constraint(m, x[8]-0.03*x[22]>=0)  #= e17: =#
    @constraint(m, x[9]-0.04*x[23]>=0)  #= e18: =#
    @constraint(m, x[10]-0.04*x[24]>=0)  #= e19: =#
    @constraint(m, x[11]-0.03*x[25]>=0)  #= e20: =#
    @constraint(m, x[12]-0.03*x[26]>=0)  #= e21: =#
    @constraint(m, x[13]-0.03*x[27]>=0)  #= e22: =#
    @constraint(m, x[14]-0.03*x[28]>=0)  #= e23: =#
    @constraint(m, x[15]-0.2*x[29]<=0)  #= e24: =#
    @constraint(m, x[16]-0.15*x[30]<=0)  #= e25: =#
    @constraint(m, x[17]<=0)  #= e26: =#
    @constraint(m, x[18]<=0)  #= e27: =#
    @constraint(m, x[19]-0.1*x[33]<=0)  #= e28: =#
    @constraint(m, x[20]-0.15*x[34]<=0)  #= e29: =#
    @constraint(m, x[21]-0.2*x[35]<=0)  #= e30: =#
    @constraint(m, x[15]-0.02*x[29]>=0)  #= e31: =#
    @constraint(m, x[16]-0.02*x[30]>=0)  #= e32: =#
    @constraint(m, x[17]-0.04*x[31]>=0)  #= e33: =#
    @constraint(m, x[18]-0.04*x[32]>=0)  #= e34: =#
    @constraint(m, x[19]-0.04*x[33]>=0)  #= e35: =#
    @constraint(m, x[20]-0.04*x[34]>=0)  #= e36: =#
    @constraint(m, x[21]-0.04*x[35]>=0)  #= e37: =#
    @constraint(m, x[22]+x[29]<=1)  #= e38: =#
    @constraint(m, x[23]+x[30]<=1)  #= e39: =#
    @constraint(m, x[24]+x[31]<=1)  #= e40: =#
    @constraint(m, x[25]+x[32]<=1)  #= e41: =#
    @constraint(m, x[26]+x[33]<=1)  #= e42: =#
    @constraint(m, x[27]+x[34]<=1)  #= e43: =#
    @constraint(m, x[28]+x[35]<=1)  #= e44: =#

    return m
end
