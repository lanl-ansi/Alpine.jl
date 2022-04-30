@testset "Utility Function Tests: Solver identifier fetch" begin

    test_solver=optimizer_with_attributes(Alpine.Optimizer, 
                          "minlp_solver" => PAVITO,
                          "nlp_solver" => IPOPT,
                          "mip_solver" => CBC,
                          "presolve_bp" => true,
                          "disc_var_pick" => 1,
                          "log_level" => 100,
                          "max_iter" => 2)

    m = blend029_gl(solver = test_solver)

    # JuMP.optimize!(m) -  # Pavito with CBC as the underlying MIP solver has issues
    alpine = JuMP.backend(m).optimizer.model

    # @test termination_status(m) == MOI.OTHER_LIMIT
    # @test alpine.logs[:n_iter] == 2

    Alpine.fetch_mip_solver_identifier(alpine;override="Pajarito.PajaritoSolver(0, Inf, 1.0e-5, false, Cbc.CbcMathProgSolverInterface.CbcSolver(Any[(:log_level, 0)]), Pajarito.UnsetSolver(), 0, false, true, Ipopt.IpoptSolver(Any[(:print_level, 0)]), true, true, false, false, true, false, false, true, false, true, true, true, true, false, true, 2.0, false, false, false, true, 1.0e-12, 1.0e-6, false, \"\")")
    @test alpine.mip_solver_id == "Pajarito"
    Alpine.fetch_mip_solver_identifier(alpine;override="Gurobi.GurobiSolver(nothing, Any[])")
    @test alpine.mip_solver_id == "Gurobi"
    Alpine.fetch_mip_solver_identifier(alpine;override="CPLEX.CplexSolver(Any[])")
    @test alpine.mip_solver_id == "Cplex"
    Alpine.fetch_mip_solver_identifier(alpine;override="Cbc.CbcMathProgSolverInterface.CbcSolver(Any[])")
    @test alpine.mip_solver_id == "Cbc"
    Alpine.fetch_mip_solver_identifier(alpine;override="GLPKMathProgInterface.GLPKInterfaceMIP.GLPKSolverMIP(false, Any[])")
    @test alpine.mip_solver_id == "GLPK"
    Alpine.fetch_mip_solver_identifier(alpine;override="Xpress.Optimizer")
    @test alpine.mip_solver_id == "Xpress"

    Alpine.fetch_nlp_solver_identifier(alpine;override="Pajarito.PajaritoSolver(0, Inf, 1.0e-5, false, Cbc.CbcMathProgSolverInterface.CbcSolver(Any[(:log_level, 0)]), Pajarito.UnsetSolver(), 0, false, true, Ipopt.IpoptSolver(Any[(:print_level, 0)]), true, true, false, false, true, false, false, true, false, true, true, true, true, false, true, 2.0, false, false, false, true, 1.0e-12, 1.0e-6, false, \"\")")
    @test alpine.nlp_solver_id == "Pajarito"
    Alpine.fetch_nlp_solver_identifier(alpine;override="Ipopt.IpoptSolver(Any[])")
    @test alpine.nlp_solver_id == "Ipopt"
    Alpine.fetch_nlp_solver_identifier(alpine;override="AmplNLWriter.AmplNLSolver(\"bonmin\", String[], \"\")")
    @test alpine.nlp_solver_id == "Bonmin"
    Alpine.fetch_nlp_solver_identifier(alpine;override="KNITRO.KnitroSolver(Any[])")
    @test alpine.nlp_solver_id == "Knitro"

    Alpine.fetch_minlp_solver_identifier(alpine;override="Pajarito.PajaritoSolver(0, Inf, 1.0e-5, false, Cbc.CbcMathProgSolverInterface.CbcSolver(Any[(:log_level, 0)]), Pajarito.UnsetSolver(), 0, false, true, Ipopt.IpoptSolver(Any[(:print_level, 0)]), true, true, false, false, true, false, false, true, false, true, true, true, true, false, true, 2.0, false, false, false, true, 1.0e-12, 1.0e-6, false, \"\")")
    @test alpine.minlp_solver_id == "Pajarito"
    Alpine.fetch_minlp_solver_identifier(alpine;override="AmplNLWriter.AmplNLSolver(\"bonmin\", String[], \"\")")
    @test alpine.minlp_solver_id == "Bonmin"
    Alpine.fetch_minlp_solver_identifier(alpine;override="KNITRO.KnitroSolver(Any[])")
    @test alpine.minlp_solver_id == "Knitro"

end

@testset "Solver Funtion Tests :: Embedding" begin
    ebdmap = Alpine.embedding_map(5)
    @test ebdmap[:L] == 2
    @test ebdmap[:H] == [[0, 0], [0, 1], [1, 1], [1, 0]]
    @test ebdmap[1] == Set([4,5])
    @test ebdmap[2] == Set([3])
    @test ebdmap[3] == Set([2,1])
    @test ebdmap[4] == Set([5,1])

    ebdmap = Alpine.embedding_map(9)
    @test ebdmap[:L] == 3
    @test ebdmap[:H_orig] == ["000", "001", "011", "010", "110", "111", "101", "100"]
    @test ebdmap[:H] == [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1], [1, 0, 0]]
    @test ebdmap[1] == Set([7,9,8,6])
    @test ebdmap[2] == Set([4,5,6])
    @test ebdmap[3] == Set([7,3])
    @test ebdmap[4] == Set([4,3,2,1])
    @test ebdmap[5] == Set([9,2,8,1])
    @test ebdmap[6] == Set([9,5,1])

    ebdmap = Alpine.embedding_map(12)
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

@testset "Utility Function Tests: check_solution_history test" begin

    test_solver=optimizer_with_attributes(Alpine.Optimizer,
    "nlp_solver" => IPOPT,
    "mip_solver" => PAVITO,
    "presolve_bt" => false,
    "disc_ratio" => 5,
    "disc_consecutive_forbid" => false)

    m1 = circle(solver = test_solver)
    JuMP.optimize!(m1)
    alpine1 = JuMP.backend(m1).optimizer.model

    Alpine.variable_values(m1)

    test_solver=optimizer_with_attributes(Alpine.Optimizer,
    "nlp_solver" => IPOPT,
    "mip_solver" => PAVITO,
    "presolve_bt" => false,
    "disc_ratio" => 5,
    "disc_consecutive_forbid" => true)
    
    m2 = circle(solver = test_solver)
    JuMP.optimize!(m2)
    alpine2 = JuMP.backend(m2).optimizer.model

    @test termination_status(m1) == MOI.OPTIMAL
    @test termination_status(m2) == MOI.OPTIMAL
    @test alpine1.logs[:n_iter] > alpine2.logs[:n_iter]
    @test isapprox(JuMP.objective_value(m1), JuMP.objective_value(m2), atol = 1E-6)

end
