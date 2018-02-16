@testset " Validation Test || Default || Integer Problem || prob03" begin

    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                       bilinear_convexhull=false,
                       monomial_convexhull=false,
                       presolve_bt=false,
                       presolve_bp=true,
                       presolve_bt_output_tol=1e-1,
                       loglevel=100)
    m = prob03(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
end
