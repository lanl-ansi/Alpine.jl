# Contains a basic model with various expressions for testing
using POD, JuMP, Ipopt, CPLEX, MathProgBase

function exprstest(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
									mip_solver=CplexSolver(),
									log_level=0))
	else
		m = Model(solver=solver)
	end

	@variable(m, px[i=1:6]>=1) # At some point if an initial value is given, keep them

	@NLconstraint(m, sum(3*px[i]^2 for i=1:4) >= 111)
	@NLconstraint(m,  -px[1] * px[2] + 4*5*px[3]*px[4] >= 222)
	@NLconstraint(m, -px[1] * px[2] <= 115)
	@NLconstraint(m, -px[1] * -px[2] >= 115)
	@NLconstraint(m, px[1] * -px[2] <= 115)
	@constraint(m, -px[1] + (-5) - 4 <= 100)
	@NLconstraint(m, px[1]+ px[2]*px[3] >= 555) # => px[1] + x23 >= 555 && x23 == px[2]*px[3]
	@NLconstraint(m, px[1]^2 - 7*px[2]^2 + px[3]^2 + px[4] <= 6666)
	@NLconstraint(m, 13*px[1] - px[2] + 5*px[3]*6 + px[4] >= 77)

	@NLobjective(m, Min, 7*px[1]*6*px[4]*2+5+17+px[1]+px[2]+px[3]+8+3*5*px[1]^2*4)

	if verbose
		print(m)
	end

	return m
end

function meyer4_expr(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=BonminNLSolver(["bonmin.iteration_limit=100"]),
									mip_solver=GurobiSolver(OutputFlag=0),
									discretization_ratio=32,
									log_level=100,
									rel_gap=0.001))
	else
		m = Model(solver=solver)
	end

    # ----- Variables ----- #
    xIdx = Int64[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119]
    @assert length(xIdx) == 119
    @variable(m, x[1:119])   # Set up all variables

    # ----> BINARY VARIABLES (IF ANY)
    binaryIdx = Int64[64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118]
    for i in binaryIdx
        setcategory(x[i], :Bin)
    end

    # ----> ALL BOUNDS
    setlowerbound(x[16], 0.0)
    setlowerbound(x[14], 0.0)
    setlowerbound(x[62], 0.0)
    setlowerbound(x[38], 0.0)
    setlowerbound(x[42], 0.0)
    setlowerbound(x[56], 0.0)
    setlowerbound(x[22], 0.0)
    setlowerbound(x[59], 0.0)
    setlowerbound(x[2], 0.0)
    setlowerbound(x[9], 0.0)
    setlowerbound(x[8], 0.0)
    setlowerbound(x[43], 0.0)
    setlowerbound(x[36], 0.0)
    setlowerbound(x[4], 0.0)
    setlowerbound(x[32], 0.0)
    setlowerbound(x[54], 0.0)
    setlowerbound(x[27], 0.0)
    setlowerbound(x[3], 0.0)
    setlowerbound(x[25], 0.0)
    setlowerbound(x[30], 0.0)
    setlowerbound(x[58], 0.0)
    setlowerbound(x[11], 0.0)
    setlowerbound(x[29], 0.0)
    setlowerbound(x[53], 0.0)
    setlowerbound(x[5], 0.0)
    setlowerbound(x[37], 0.0)
    setlowerbound(x[63], 0.0)
    setlowerbound(x[57], 0.0)
    setlowerbound(x[55], 0.0)
    setlowerbound(x[24], 0.0)
    setlowerbound(x[41], 0.0)
    setlowerbound(x[18], 0.0)
    setlowerbound(x[52], 0.0)
    setlowerbound(x[1], 0.0)
    setlowerbound(x[7], 0.0)
    setlowerbound(x[13], 0.0)
    setlowerbound(x[49], 0.0)
    setlowerbound(x[21], 0.0)
    setlowerbound(x[10], 0.0)
    setlowerbound(x[26], 0.0)
    setlowerbound(x[45], 0.0)
    setlowerbound(x[12], 0.0)
    setlowerbound(x[40], 0.0)
    setlowerbound(x[44], 0.0)
    setlowerbound(x[61], 0.0)
    setlowerbound(x[50], 0.0)
    setlowerbound(x[31], 0.0)
    setlowerbound(x[33], 0.0)
    setlowerbound(x[47], 0.0)
    setlowerbound(x[28], 0.0)
    setlowerbound(x[35], 0.0)
    setlowerbound(x[6], 0.0)
    setlowerbound(x[60], 0.0)
    setlowerbound(x[17], 0.0)
    setlowerbound(x[23], 0.0)
    setlowerbound(x[34], 0.0)
    setlowerbound(x[46], 0.0)
    setlowerbound(x[51], 0.0)
    setlowerbound(x[19], 0.0)
    setlowerbound(x[48], 0.0)
    setlowerbound(x[20], 0.0)
    setlowerbound(x[39], 0.0)
    setlowerbound(x[15], 0.0)
    setupperbound(x[1],20.0)
    setupperbound(x[2],20.0)
    setupperbound(x[3],20.0)
    setupperbound(x[4],20.0)
    setupperbound(x[5],50.0)
    setupperbound(x[6],50.0)
    setupperbound(x[7],50.0)
    setupperbound(x[8],50.0)
    setupperbound(x[9],47.5)
    setupperbound(x[10],47.5)
    setupperbound(x[11],47.5)
    setupperbound(x[12],47.5)
    setupperbound(x[13],28.0)
    setupperbound(x[14],28.0)
    setupperbound(x[15],28.0)
    setupperbound(x[16],28.0)
    setupperbound(x[17],100.0)
    setupperbound(x[18],100.0)
    setupperbound(x[19],100.0)
    setupperbound(x[20],100.0)
    setupperbound(x[21],30.0)
    setupperbound(x[22],30.0)
    setupperbound(x[23],30.0)
    setupperbound(x[24],30.0)
    setupperbound(x[25],25.0)
    setupperbound(x[26],25.0)
    setupperbound(x[27],25.0)
    setupperbound(x[28],25.0)
    setupperbound(x[29],300.5)
    setupperbound(x[30],300.5)
    setupperbound(x[31],300.5)
    setupperbound(x[32],300.5)
    setupperbound(x[33],300.5)
    setupperbound(x[34],300.5)
    setupperbound(x[35],300.5)
    setupperbound(x[36],300.5)
    setupperbound(x[37],300.5)
    setupperbound(x[38],300.5)
    setupperbound(x[39],300.5)
    setupperbound(x[40],300.5)
    setupperbound(x[41],300.5)
    setupperbound(x[42],300.5)
    setupperbound(x[43],300.5)
    setupperbound(x[44],300.5)
    setupperbound(x[45],20.0)
    setupperbound(x[46],50.0)
    setupperbound(x[47],47.5)
    setupperbound(x[48],28.0)
    setupperbound(x[49],100.0)
    setupperbound(x[50],30.0)
    setupperbound(x[51],25.0)
    setupperbound(x[52],12.0)
    setupperbound(x[53],175.0)
    setupperbound(x[54],100.0)
    setupperbound(x[55],1200.0)
    setupperbound(x[56],227.5)
    setupperbound(x[57],200.0)
    setupperbound(x[58],1080.0)
    setupperbound(x[59],17.5)
    setupperbound(x[60],2000.0)
    setupperbound(x[61],360.0)
    setupperbound(x[62],1400.0)
    setupperbound(x[63],1400.0)

    # ----> OBJECTIVE FUNCTION
    @objective(m, Min, x[119])

    # ----> NON-LINEAR CONSTRAINTS
    @NLconstraint(m, e26,0.01*(x[55]*x[36]+x[58]*x[39]+x[61]*x[42])-(x[52]*x[29]+x[52]*x[33]+x[52]*x[34]+x[52]*x[35])+x[1]+8.00000000000001*x[5]+4*x[9]+12*x[13]+5*x[17]+0.5*x[21]+10*x[25]==0.0)
    @NLconstraint(m, e27,800*(0.1*(x[56]*x[36]+x[59]*x[39]+x[62]*x[42])-(x[53]*x[29]+x[53]*x[33]+x[53]*x[34]+x[53]*x[35]))+50*x[1]+175*x[5]+8*x[9]+100*x[13]+70*x[17]+10*x[21]+5*x[25]==0.0)
    @NLconstraint(m, e28,0.05*(x[57]*x[36]+x[60]*x[39]+x[63]*x[42])-0.02*(x[54]*x[29]+x[54]*x[33]+x[54]*x[34]+x[54]*x[35])+25*x[1]+100*x[5]+5*x[9]+20*x[13]+12.5*x[17]+2.5*x[21]+7.50000000000001*x[25]==0.0)
    @NLconstraint(m, e38,x[52]*x[29]+x[55]*x[30]+x[58]*x[31]+x[61]*x[32]-5*x[29]-5*x[30]-5*x[31]-5*x[32]+95*x[45]+795*x[46]+395*x[47]+1195*x[48]+495*x[49]+45*x[50]+995*x[51]<=0.0)

    @constraint(m, e1,-75.0708333333333*x[1]-150.141666666667*x[2]-280.264444444444*x[3]-245.231388888889*x[4]-55.0519444444444*x[5]-125.118055555556*x[6]-260.245555555556*x[7]-215.203055555556*x[8]-30.0283333333333*x[9]-115.108611111111*x[10]-240.226666666667*x[11]-220.207777777778*x[12]-55.0519444444444*x[13]-140.132222222222*x[14]-245.231388888889*x[15]-245.231388888889*x[16]-55.0519444444444*x[17]-40.0377777777778*x[18]-150.141666666667*x[19]-150.141666666667*x[20]-40.0377777777778*x[21]-120.113333333333*x[22]-230.217222222222*x[23]-230.217222222222*x[24]-30.0283333333333*x[25]-60.0566666666667*x[26]-175.165277777778*x[27]-165.155833333333*x[28]-1177.97083333333*x[29]-2975.27555555556*x[30]-1263.05111111111*x[31]-1293.07944444444*x[32]-1182.97555555556*x[33]-1313.09833333333*x[34]-1293.07944444444*x[35]-2975.27555555556*x[36]-3025.32277777778*x[37]-2995.29444444444*x[38]-1313.09833333333*x[39]-1233.02277777778*x[40]-1213.00388888889*x[41]-1293.07944444444*x[42]-1202.99444444444*x[43]-1213.00388888889*x[44]-150.141666666667*x[45]-135.1275*x[46]-100.094444444444*x[47]-90.085*x[48]-40.0377777777778*x[49]-70.0661111111111*x[50]-45.0425*x[51]-9345*x[64]-18690*x[65]-34888*x[66]-30527*x[67]-6853*x[68]-15575*x[69]-32396*x[70]-26789*x[71]-3738*x[72]-14329*x[73]-29904*x[74]-27412*x[75]-6853*x[76]-17444*x[77]-30527*x[78]-30527*x[79]-6853*x[80]-4984*x[81]-18690*x[82]-18690*x[83]-4984*x[84]-14952*x[85]-28658*x[86]-28658*x[87]-3738*x[88]-7476*x[89]-21805*x[90]-20559*x[91]-9345*x[92]-9968*x[93]-19936*x[94]-23674*x[95]-9968*x[96]-26166*x[97]-23674*x[98]-9968*x[99]-16198*x[100]-12460*x[101]-26166*x[102]-16198*x[103]-13706*x[104]-23674*x[105]-12460*x[106]-13706*x[107]-18690*x[108]-16821*x[109]-12460*x[110]-11214*x[111]-4984*x[112]-8722*x[113]-5607*x[114]-13972*x[115]-36676*x[116]-13972*x[117]-13972*x[118]+x[119]==0.0)
    @constraint(m, e2,-x[1]-x[2]-x[3]-x[4]-x[45]<=-20.0)
    @constraint(m, e140,x[35]-300.5*x[98]<=0.0)

    if verbose
        print(m)
    end

    return m
end

function operator_b(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(OutputFlag=0),
								   presolve_bound_tightening=true,
								   presolve_bound_tightening_algo=2,
								   presolve_bt_output_tol=1e-1,
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	@variable(m, x[1:4]>=0)
	@variable(m, y[1:3]<=0)
																# PARSING TARGETS #
	@constraint(m, x[1] + x[2] + -y[1] >= 1)					# None
	@constraint(m, 5x[1] + 8x[2] + 9y[2] <= 2)					# None
	@constraint(m, 4*x[1] + 5*x[2] + -3*y[1] >= 3)				# None
	@constraint(m, 3 <= x[1] + x[2] <= 4)						# None

	@NLconstraint(m, x[1]*x[2] + 3*x[2]*x[3] >= 5)				# x[1]*x[2] and 3*x[2]*x[3]
	@NLconstraint(m, x[1]*(x[2]+4) + 10*44*x[2]*5*x[3] <= 6)	# x[2]*x[3]
	@NLconstraint(m, x[1]*4 + x[2]*(x[2]+x[3]) >= 7)			# None
	@NLconstraint(m, 3*x[1]*x[2] + x[2]*-x[3] <= 8)				# x[1]*x[2]

	@NLconstraint(m, -1*x[1]^2 - x[2]^2 >= 9)					# x[1]^2 and x[2]^2
	@NLconstraint(m, x[1]^(1+1) + 5*x[2]^2 + 3*x[3]^3 <= 10)	# x[1]^2 and x[2]^2
	@NLconstraint(m, (x[1]+5)^2 >= 11)							# None
	@NLconstraint(m, x[1]^3 + x[3]^99 >= 12)					# None

	@NLconstraint(m, x[1]*x[2]*x[3] <= 13)						# x[1]*x[2]*x[3]
	@NLconstraint(m, (x[1]*x[2])*x[3] >= 14)					# x[1]*x[2]  ******
	@NLconstraint(m, x[1]*(x[2]*x[3]) <= 15)					# x[2]*x[3]  ******
	@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 16)					# x[1]*x[2]*x[3]*x[4]
	@NLconstraint(m, (x[1]*x[2])*(x[3]*x[4]) <= 17)				# x[1]*x[2] and x[3]*x[4]
	@NLconstraint(m, x[1]*(x[2]*x[3])*x[4] >= 18)				# x[2]*x[3]  ******
	@NLconstraint(m, (x[1]*x[2])*x[3]*x[4] <= 19)				# x[1]*x[2]  ******
	@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 20)				# x[1]*x[2]	 ******
	@NLconstraint(m, (x[1]*(x[2]*x[3])*x[4]) <= 21)				# x[2]*x[3]	 ******

	@NLconstraint(m, 4*5*6*x[1]*x[2]*x[3]*x[4] >= 22)			# x[1]*x[2]*x[3]*x[4]
	@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4] <= 23)			# None
	@NLconstraint(m, x[1]*x[1]*x[2]*x[2]*x[3]*x[3] >= 24)		# x[1]*x[1]*x[2]*x[2]*x[3]*x[3]
	@NLconstraint(m, (x[1]+1)*(x[2]+2)*(x[3]+3)*(x[4]+4) <= 25)	# None

	@NLconstraint(m, 50sin(x[1]) - 32*cos(x[2]) >= 26)			# sin(x[1]), cos(x[2])
	@NLconstraint(m, sin(x[1]+x[2]) - cos(x[2]-x[3]) + sin(-x[2]+x[3]) <= 27)
	@NLconstraint(m, sin(4*x[1]*x[2]) + cos(x[2]*-1*x[3]) >= 28)# x[1]*x[2] and x[2]*-1*x[3]  ******
	@NLconstraint(m, sin(x[1]*x[2]*x[3]) <= 29)					# x[1]*x[2]*x[3]

	if verbose
		print(m)
	end

	return m
end

function operator_basic(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(),
								   log_level=1))
	else
		m = Model(solver=solver)
	end

	@variable(m, x[1:4]>=0)

	@NLconstraint(m, x[1]^2 >= 1)
	@NLconstraint(m, x[1]*x[2] <= 1)
	@NLconstraint(m, x[1]^2 + x[2]*x[3] <= 1)

	@NLconstraint(m, x[1] * (x[2]*x[3]) >= 1)
	@NLconstraint(m, x[1]^2 * (x[2]^2 * x[3]^2) <= 1)

	@NLconstraint(m, (x[1] * x[2]) * x[3] >= 1)
	@NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

	@NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
	@NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2 <= 1)
	@NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
	@NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

	@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 1)
	@NLconstraint(m, ((x[1]^2*x[2])*x[3])*x[4] <= 1)
	@NLconstraint(m, ((x[1]*x[2]^2)*x[3])*x[4] >= 1)
	@NLconstraint(m, ((x[1]*x[2])*x[3]^2)*x[4] <= 1)
	@NLconstraint(m, ((x[1]*x[2])*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, ((x[1]^2*x[2]^2)*x[3]^2)*x[4]^2 <=1)

	@NLconstraint(m, x[1]*(x[2]*(x[3]*x[4])) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*(x[3]*x[4])) <= 1)
	@NLconstraint(m, x[1]*(x[2]^2*(x[3]*x[4])) >= 1)
	@NLconstraint(m, x[1]*(x[2]*(x[3]^2*x[4])) <= 1)
	@NLconstraint(m, x[1]*(x[2]*(x[3]*x[4]^2)) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]^2*(x[3]^2*x[4]^2)) <= 1)

	@NLconstraint(m, x[1]*x[2]*x[3] >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3] >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3] >= 1)
	@NLconstraint(m, x[1]*x[2]*x[3]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*x[3] >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2 >= 1)

	@NLconstraint(m, (x[1]*x[2])*(x[3]*x[4]) >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]) >= 1)
	@NLconstraint(m, (x[1]*x[2]^2)*(x[3]*x[4]) >= 1)
	@NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]) >= 1)
	@NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]^2) >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]^2) >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*(x[3]^2*x[4]) >= 1)

	@NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3]*x[4] >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3]*x[4] >= 1)
	@NLconstraint(m, x[1]*x[2]*x[3]^2*x[4]^2 >= 1)
	@NLconstraint(m, x[1]*x[2]^2*x[3]^2*x[4] >= 1)
	@NLconstraint(m, x[1]^2*x[2]*x[3]*x[4]^2 >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4]^2 >= 1)

	@NLconstraint(m, (x[1]*x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, (x[1]^2*x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]^2*x[3])*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]*x[3]^2)*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, (x[1]^2*x[2]^2*x[3]^2)*x[4]^2 >= 1)

	@NLconstraint(m, x[1]*(x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]^2*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4]^2 >= 1)
	@NLconstraint(m, x[1]*(x[2]^2*x[3]^2)*x[4] >= 1)
	@NLconstraint(m, x[1]^2*(x[2]^2*x[3])*x[4] >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4]^2 >= 1)

	@NLconstraint(m, x[1]*(x[2]*x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*(x[2]^2*x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]^2*x[4]) >= 1)
	@NLconstraint(m, x[1]*(x[2]*x[3]*x[4]^2) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]*x[3]^2*x[4]) >= 1)
	@NLconstraint(m, x[1]^2*(x[2]^2*x[3]^2*x[4]^2) >= 1)

	@NLconstraint(m, x[1]*x[2]*(x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]^2*x[2]*(x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*x[2]^2*(x[3]*x[4]) >= 1)
	@NLconstraint(m, x[1]*x[2]*(x[3]^2*x[4]) >= 1)
	@NLconstraint(m, x[1]*x[2]*(x[3]*x[4]^2) >= 1)
	@NLconstraint(m, x[1]^2*x[2]^2*(x[3]^2*x[4]^2) >= 1)

	@NLconstraint(m, (x[1]*x[2])*x[3]*x[4] >= 1)
	@NLconstraint(m, (x[1]^2*x[2])*x[3]*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2]^2)*x[3]*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4] >= 1)
	@NLconstraint(m, (x[1]*x[2])*x[3]*x[4]^2 >= 1)
	@NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4]^2 >= 1)
	@NLconstraint(m, (x[1]^2*x[2]^2)*x[3]^2*x[4]^2 >= 1)

	if verbose
		print(m)
	end

	return m
end

function operator_c(;verbose=false, solver=nothing)

	if solver == nothing
		m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
							   	   mip_solver=CplexSolver()))
	else
		m = Model(solver=solver)
	end

	@variable(m, px[i=1:6]>=1) # At some point if an initial value is given, keep them

	@NLconstraint(m, sum(3*px[i]^2 for i=1:4) >= 111)
	@NLconstraint(m,  -1 * px[1]*px[2] + 4*5*px[3]*px[4] >= 222)
	@NLconstraint(m, px[1]+ px[2] * px[3] >= 555) # => px[1] + x23 >= 555 && x23 == px[2]*px[3]
	@NLconstraint(m, px[1]^2 - 7 * px[2]^2 + px[3]^2 + px[4] <= 6666)
	@NLconstraint(m, 13*px[1] - px[2] + 5*px[3]*6 + px[4] >= 77)

	@NLobjective(m, Min, 7*px[1]*6*px[4]*2+5+17+px[1]+px[2]+px[3]+8+3*5*px[1]^2*4)
	if verbose
		print(m)
	end

	return m
end
