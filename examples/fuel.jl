using POD, JuMP, CPLEX, Ipopt, MathProgBase, AmplNLWriter, CoinOptServices

function castro02m2(;verbose=false, solver=nothing)

    if solver == nothing
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    								mip_solver=CplexSolver(CPX_PARAM_SCRIND=0),
                                    bilinear_convexhull=false,
                                    presolve_bound_tightening=false,
                                    log_level=100,
                                    rel_gap=0.001))
    else
        m = Model(solver=solver)
    end

    #define variables and bounds
    @variable(m, x[1:15])

    setcategory(x[1], :Bin)
    setcategory(x[2], :Bin)
    setcategory(x[3], :Bin)

    # Bounds Constraints
    @constraint(m, x[7] >= 50.0)
    @constraint(m, x[7] <= 700.0)
    @constraint(m, x[8] >= 50.0)
    @constraint(m, x[8] <= 700.0)
    @constraint(m, x[9] >= 50.0)
    @constraint(m, x[9] <= 700.0)
    @constraint(m, x[10] >= 0.0)
    @constraint(m, x[11] >= 0.0)
    @constraint(m, x[12] >= 0.0)
    @constraint(m, x[13] >= 0.0)
    @constraint(m, x[13] <= 4000.0)
    @constraint(m, x[14] >= 0.0)
    @constraint(m, x[14] <= 4000.0)
    @constraint(m, x[15] >= 2000.0)
    @constraint(m, x[15] <= 4000.0)

    # linear constraints
    @constraint(m, -100*x[1]+x[4]>=0)  #= e2: =#
    @constraint(m, -100*x[2]+x[5]>=0)  #= e3: =#
    @constraint(m, -100*x[3]+x[6]>=0)  #= e4: =#
    @constraint(m, -500*x[1]+x[4]<=0)  #= e5: =#
    @constraint(m, -500*x[2]+x[5]<=0)  #= e6: =#
    @constraint(m, -500*x[3]+x[6]<=0)  #= e7: =#
    @constraint(m, x[10]+x[13] ==3500)  #= e8: =#
    @constraint(m, x[11]-x[13]+x[14] ==500)  #= e9: =#
    @constraint(m, x[12]-x[14]+x[15] ==500)  #= e10: =#
    @constraint(m, x[4]+x[7]>=400)  #= e11: =#
    @constraint(m, x[5]+x[8]>=900)  #= e12: =#
    @constraint(m, x[6]+x[9]>=700)  #= e13: =#

    @objective(m, Min, 0.0025*x[7]^2 + 6*x[7] + 0.0025*x[8]^2 + 6*x[8] + 0.0025*x[9]^2 + 6*x[9] + 900)
    @NLconstraint(m, -(0.005*x[4]^2+x[4])-50*x[1]+x[10] ==0)  #= e14: =#
    @NLconstraint(m, -(0.005*x[5]^2+x[5])-50*x[2]+x[11] ==0)  #= e15: =#
    @NLconstraint(m, -(0.005*x[6]^2+x[6])-50*x[3]+x[12] ==0)  #= e16: =#

    if verbose
        print(m)
    end

    return m
end
