@testset "Utility Function Tests: Solver identifier fetch" begin

    test_solver=PODSolver(minlp_local_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                          nlp_local_solver=IpoptSolver(print_level=0),
                          mip_solver=CbcSolver(logLevel=0),
                          presolve_bp=true,
                          disc_var_pick=1,
                          log=100,
                          maxiter=2,
                          presolve_bt_width_tol=1e-3,
                          presolve_bt=false)
    m = blend029_gl(solver=test_solver)
    status = solve(m)

    @test status == :Infeasible
    @test m.internalModel.logs[:n_iter] == 2

    POD.fetch_mip_solver_identifier(m.internalModel;override="Pajarito.PajaritoSolver(0, Inf, 1.0e-5, false, Cbc.CbcMathProgSolverInterface.CbcSolver(Any[(:logLevel, 0)]), Pajarito.UnsetSolver(), 0, false, true, Ipopt.IpoptSolver(Any[(:print_level, 0)]), true, true, false, false, true, false, false, true, false, true, true, true, true, false, true, 2.0, false, false, false, true, 1.0e-12, 1.0e-6, false, \"\")")
    @test m.internalModel.mip_solver_identifier == "Pajarito"
    POD.fetch_mip_solver_identifier(m.internalModel;override="Gurobi.GurobiSolver(nothing, Any[])")
    @test m.internalModel.mip_solver_identifier == "Gurobi"
    POD.fetch_mip_solver_identifier(m.internalModel;override="CPLEX.CplexSolver(Any[])")
    @test m.internalModel.mip_solver_identifier == "CPLEX"
    POD.fetch_mip_solver_identifier(m.internalModel;override="Cbc.CbcMathProgSolverInterface.CbcSolver(Any[])")
    @test m.internalModel.mip_solver_identifier == "Cbc"
    POD.fetch_mip_solver_identifier(m.internalModel;override="GLPKMathProgInterface.GLPKInterfaceMIP.GLPKSolverMIP(false, Any[])")
    @test m.internalModel.mip_solver_identifier == "GLPK"

    POD.fetch_nlp_solver_identifier(m.internalModel;override="Pajarito.PajaritoSolver(0, Inf, 1.0e-5, false, Cbc.CbcMathProgSolverInterface.CbcSolver(Any[(:logLevel, 0)]), Pajarito.UnsetSolver(), 0, false, true, Ipopt.IpoptSolver(Any[(:print_level, 0)]), true, true, false, false, true, false, false, true, false, true, true, true, true, false, true, 2.0, false, false, false, true, 1.0e-12, 1.0e-6, false, \"\")")
    @test m.internalModel.nlp_local_solver_identifier == "Pajarito"
    POD.fetch_nlp_solver_identifier(m.internalModel;override="Ipopt.IpoptSolver(Any[])")
    @test m.internalModel.nlp_local_solver_identifier == "Ipopt"
    POD.fetch_nlp_solver_identifier(m.internalModel;override="AmplNLWriter.AmplNLSolver(\"bonmin\", String[], \"\")")
    @test m.internalModel.nlp_local_solver_identifier == "Bonmin"
    @test "NLopt" == "NLopt"
    POD.fetch_nlp_solver_identifier(m.internalModel;override="KNITRO.KnitroSolver(Any[])")
    @test m.internalModel.nlp_local_solver_identifier == "Knitro"

    POD.fetch_minlp_solver_identifier(m.internalModel;override="Pajarito.PajaritoSolver(0, Inf, 1.0e-5, false, Cbc.CbcMathProgSolverInterface.CbcSolver(Any[(:logLevel, 0)]), Pajarito.UnsetSolver(), 0, false, true, Ipopt.IpoptSolver(Any[(:print_level, 0)]), true, true, false, false, true, false, false, true, false, true, true, true, true, false, true, 2.0, false, false, false, true, 1.0e-12, 1.0e-6, false, \"\")")
    @test m.internalModel.minlp_local_solver_identifier == "Pajarito"
    POD.fetch_minlp_solver_identifier(m.internalModel;override="AmplNLWriter.AmplNLSolver(\"bonmin\", String[], \"\")")
    @test m.internalModel.minlp_local_solver_identifier == "Bonmin"
    POD.fetch_minlp_solver_identifier(m.internalModel;override="KNITRO.KnitroSolver(Any[])")
    @test m.internalModel.minlp_local_solver_identifier == "Knitro"
    @test "NLopt" == "NLopt"
end

@testset "Solver Funtion Tests :: Embedding" begin
    ebdmap = POD.embedding_map(5)
    @test ebdmap[:L] == 2
    @test ebdmap[:H] == [[0, 0], [0, 1], [1, 1], [1, 0]]
    @test ebdmap[1] == Set([4,5])
    @test ebdmap[2] == Set([3])
    @test ebdmap[3] == Set([2,1])
    @test ebdmap[4] == Set([5,1])

    ebdmap = POD.embedding_map(9)
    @test ebdmap[:L] == 3
    @test ebdmap[:H_orig] == ["000", "001", "011", "010", "110", "111", "101", "100"]
    @test ebdmap[:H] == [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1], [1, 0, 0]]
    @test ebdmap[1] == Set([7,9,8,6])
    @test ebdmap[2] == Set([4,5,6])
    @test ebdmap[3] == Set([7,3])
    @test ebdmap[4] == Set([4,3,2,1])
    @test ebdmap[5] == Set([9,2,8,1])
    @test ebdmap[6] == Set([9,5,1])

    ebdmap = POD.embedding_map(12)
    @test ebdmap[:L] == 4
    @test ebdmap[:H_orig] == ["0000", "0001", "0011", "0010", "0110", "0111", "0101", "0100", "1100", "1101", "1111", "1110", "1010", "1011", "1001", "1000"]
    @test ebdmap[:H] == [[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 1, 0], [0, 1, 1, 0], [0, 1, 1, 1], [0, 1, 0, 1], [0, 1, 0, 0], [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 1, 1], [1, 1, 1, 0], [1, 0, 1, 0], [1, 0, 1, 1], [1, 0, 0, 1], [1, 0, 0, 0]]
    @test ebdmap[1] == Set([10,11,12])
    @test ebdmap[2] == Set([7,9,10,11,8,6,12])
    @test ebdmap[3] == Set([4,5,6,12])
    @test ebdmap[4] == Set([7,3,11,12])
    @test ebdmap[5] == Set([7,4,2,3,5,8,6,1])
    @test ebdmap[6] == Set([4,2,3,1])
    @test ebdmap[7] == Set([9,10,2,8,1])
    @test ebdmap[8] == Set([9,5,1])
end
