@testset " Validation Test || Default || Integer || prob03" begin

    test_solver = PODSolver(minlp_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                            nlp_solver=IpoptSolver(print_level=0),
                            mip_solver=CbcSolver(logLevel=0),
                            int2bin=true,
                            loglevel=100)

    m = st_miqp4(solver=test_solver)
    status = solve(m)
    @test status == :Optimal
    @test isapprox(m.objVal, -4574.0; atol=1e-4)
    @test isapprox(m.objBound, -4574.0; atol=1e-4)
    @test m.colVal[1] == -4574.0
    @test m.colVal[2] == 5.0
    @test m.colVal[3] == 10.0
    @test m.colVal[4] == 15.0
    @test m.colVal[5] == 1.0
    @test m.colVal[6] == 1.0
    @test m.colVal[7] == 1.0
    @test m.internalModel.logs[:n_iter] == 1
end
