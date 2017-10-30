@testset "Solving algorithm tests" begin

    @testset " Validation Test || AMP-TMC || basic solve || exampls/nlp1.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                           bilinear_convexhull=false,
                           monomial_convexhull=false,
                           presolve_bound_tightening=false,
                           presolve_bt_output_tol=1e-1,
                           log_level=0)
        m = nlp1(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 7

        # Out-dated tests base on initial implementation [PASSED on previouse tests]
        # @test isapprox(m.internalModel.logs[:obj][1], 58.38367118523376; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:obj][2], 58.38367118523376; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:obj][3], 58.38367118523376; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:obj][4], 58.38367118523376; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:obj][5], 58.38367118523376; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:obj][6], 58.38367118523376; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:obj][7], 58.38367118523376; atol=1e-3)
        #
        # @test isapprox(m.internalModel.logs[:bound][1], 25.7091; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][2], 51.5131; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][3], 56.5857; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][4], 58.2394; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][5], 58.3611; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][6], 58.3681; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][7], 58.3824; atol=1e-3)

        # Out-dated tests base on initial implementation [PASSED]
        # Testing discretizations for x1 and x2
        # ans = Dict()
        # ans[1] = [1.0,1.12035,2.10076,2.43359,2.51991,2.52774,2.55118,2.60884,2.64305,2.75053,2.89483,3.02324,4.80577, 10.0]
        # ans[2] = [1.0,1.13852,2.26689,2.64995,2.89107,3.02633,3.06707,3.13215,3.15649,3.18081,3.3286,5.38017, 10.0]
        #
        # @test length(ans[1]) == length(m.internalModel.discretization[1])
        # @test length(ans[2]) == length(m.internalModel.discretization[2])
        # for i in 1:length(ans[1])
        #     @test isapprox(ans[1][i], m.internalModel.discretization[1][i]; atol=1e-1)
        # end
        # for i in 1:length(ans[2])
        #     @test isapprox(ans[2][i], m.internalModel.discretization[2][i]; atol=1e-1)
        # end
    end

    @testset " Validation Test || AMP-TMC || basic solve || examples/nlp3.jl (3 iterations)" begin

        test_solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
        					   mip_solver=CbcSolver(logLevel=0),
                               bilinear_convexhull=false,
                               monomial_convexhull=false,
        					   log_level=0,
                               maxiter=3,
        					   presolve_bt_width_tol=1e-3,
        					   presolve_bound_tightening=false,
        					   discretization_var_pick_algo=0)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits

        @test isapprox(m.objVal, 7049.247897696512; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 3

        # Out-dated tests base on initial implementation [PASSED]
        # @test isapprox(m.internalModel.logs[:bound][1], 3183.6748; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][2], 5105.0043; atol=1e-3)
        # @test isapprox(m.internalModel.logs[:bound][3], 6485.8829; atol=1e-3)
        #
        # ans = Dict()
        # ans[1] = [100.0, 424.721, 838.577, 3054.31, 10000.0]
        # ans[2] = [1000.0, 1216.0, 1864.0, 3609.97, 10000.0]
        # ans[3] = [1000.0, 1407.17, 2337.16, 2859.97, 4989.93, 7359.97, 10000.0]
        # ans[4] = [10.0, 87.0177, 132.916, 168.929, 277.018, 390.0]
        # ans[5] = [10.0, 103.101, 185.828, 214.871, 311.121, 378.328, 488.101, 780.0]
        # ans[6] = [10.0, 122.982, 231.071, 267.084, 312.982, 390.0]
        # ans[7] = [10.0, 93.9165, 143.1, 242.272, 335.6, 478.917, 780.0]
        # ans[8] = [10.0, 178.101, 273.328, 308.621, 417.371, 490.828, 613.101, 880.0]
        #
        # for i in 1:8
        #     for j in 1:length(ans[i])
        #         @test length(ans[i]) == length(m.internalModel.discretization[i])
        #         @test isapprox(ans[i][j], m.internalModel.discretization[i][j]; atol=1e-1)
        #     end
        # end
    end

    # Diabled due to solver scheme mismatch (PajaritoSolver doesn't take pure continuous problem)
    # @testset " Validation Test || BT-AMP-TMC || basic solve || examples/nlp1.jl " begin
    #
    #     # PASSED test but disabled due to no open-source solver availibility
    #
    #     test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    #                            mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(), log_level=0),
    #                            bilinear_convexhull=false,
    #                            monomial_convexhull=false,
    #                            presolve_bound_tightening=true,
    #                            presolve_bound_tightening_algo=1,
    #                            presolve_bt_output_tol=1e-1,
    #                            log_level=0)
    #
    #     m = nlp1(solver=test_solver)
    #     @variable(m, red_var, Bin)
    #     status = solve(m)
    #
    #     @test status == :Optimal
    #     @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    #     @test m.internalModel.logs[:n_iter] == 3
    #
    #     # Out-dated tests base on initial implementation [PASSED]
    #     # @test isapprox(m.internalModel.l_var_tight[1], 2.4; atol=1e-2)
    #     # @test isapprox(m.internalModel.l_var_tight[2], 3.0; atol=1e-2)
    #     # @test isapprox(m.internalModel.l_var_tight[3], 5.76; atol=1e-2)
    #     # @test isapprox(m.internalModel.l_var_tight[4], 9.0; atol=1e-2)
    #     # @test isapprox(m.internalModel.l_var_tight[5], 8.0; atol=1e-2)
    #     #
    #     # @test isapprox(m.internalModel.u_var_tight[1], 2.7; atol=1e-2)
    #     # @test isapprox(m.internalModel.u_var_tight[2], 3.3; atol=1e-2)
    #     # @test isapprox(m.internalModel.u_var_tight[3], 7.29; atol=1e-2)
    #     # @test isapprox(m.internalModel.u_var_tight[4], 10.89; atol=1e-2)
    #     # @test isapprox(m.internalModel.u_var_tight[5], 8.91; atol=1e-2)
    # end

    @testset " Validation Test || PBT-AMP-TMC || basic solve || exampls/nlp1.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    							   mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                                   bilinear_convexhull=false,
                                   monomial_convexhull=false,
    							   presolve_bound_tightening=true,
    							   presolve_bound_tightening_algo=2,
                                   presolve_bt_output_tol=1e-1,
    							   log_level=0)
        m = nlp1(solver=test_solver)
        status = solve(m)

        # @show m.internalModel.l_var_tight
        # @show m.internalModel.u_var_tight

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 2

        # Out-dated tests base on initial implementation [PASSED]
        # @show m.internalModel.l_var_tight
        # @test isapprox(m.internalModel.l_var_tight[1], 2.5; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[2], 3.1; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[3], 6.25; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[4], 9.61; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[5], 8.0; atol=1e-2)
        #
        # @test isapprox(m.internalModel.u_var_tight[1], 2.6; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[2], 3.2; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[3], 6.76; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[4], 10.24; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[5], 8.32; atol=1e-2)
    end

    @testset " Validation Test || BT-AMP-TMC || basic solve || examples/nlp3.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    							   mip_solver=CbcSolver(logLevel=0),
                                   bilinear_convexhull=false,
    							   log_level=0,
                                   maxiter=3,
    							   presolve_bt_width_tol=1e-3,
    							   presolve_bt_output_tol=1e-1,
    							   presolve_bound_tightening=true,
                                   presolve_bound_tightening_algo=1,
    							   presolve_maxiter=2,
    							   discretization_var_pick_algo=max_cover_var_picker)
        m = nlp3(solver=test_solver)

        status = solve(m)

        @test status == :UserLimits
        @test m.internalModel.logs[:n_iter] == 3
        @test m.internalModel.logs[:bt_iter] == 2

        # @test isapprox(m.internalModel.u_var_tight[1], 4640.4; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[2], 5634.2; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[3], 5920.1; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[4], 334.2; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[5], 588.4; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[6], 390.0; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[7], 626.3; atol=1e-2)
        # @test isapprox(m.internalModel.u_var_tight[8], 688.4; atol=1e-2)
        #
        # @test isapprox(m.internalModel.l_var_tight[1], 100.0; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[2], 1000.0; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[3], 1000.0; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[4], 10.0; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[5], 107.8; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[6], 10.0; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[7], 16.9; atol=1e-2)
        # @test isapprox(m.internalModel.l_var_tight[8], 89.2; atol=1e-2)
    end

    @testset " Validation Test || PBT-AMP-TMC || basic solve || examples/nlp3.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                   mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(), log_level=0),
                                   bilinear_convexhull=false,
                                   log_level=0,
                                   maxiter=2,
                                   presolve_bound_tightening=true,
                                   presolve_bt_width_tol=1e-3,
                                   presolve_bt_output_tol=1e-1,
                                   presolve_bound_tightening_algo=2,
                                   presolve_maxiter=2,
                                   discretization_var_pick_algo=max_cover_var_picker)

        m = nlp3(solver=test_solver)
        status = solve(m)
        @test status == :UserLimits
        @test m.internalModel.logs[:n_iter] == 2
        @test m.internalModel.logs[:bt_iter] == 2

        # @test isapprox(m.internalModel.u_var_tight[1], 3011.9; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[2], 3951.3; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[3], 5836.4; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[4], 262.1; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[5], 390.5; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[6], 388.0; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[7], 370.9; atol=1e-1)
        # @test isapprox(m.internalModel.u_var_tight[8], 490.5; atol=1e-1)
        #
        # @test isapprox(m.internalModel.l_var_tight[1], 100.0; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[2], 1000.0; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[3], 2241.3; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[4], 12.0; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[5], 252.0; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[6], 10.0; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[7], 100.2; atol=1e-1)
        # @test isapprox(m.internalModel.l_var_tight[8], 352.0; atol=1e-1)
    end

    @testset " Validation Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(), log_level=0),
                           bilinear_convexhull=true,
                           monomial_convexhull=true,
                           presolve_bound_tightening=false,
                           log_level=0)
        m = nlp1(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 7
    end

    # Diabled due to solver scheme mismatch (PajaritoSolver doesn't take pure continuous problem)
    # @testset " Validation Test || BT-AMP-CONV || basic solve || examples/nlp1.jl" begin
    #     test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    #                        mip_solver=GurobiSolver(OutputFlag=0),
    #                        bilinear_convexhull=true,
    #                        monomial_convexhull=true,
    #                        presolve_bound_tightening=true,
    #                        presolve_bound_tightening_algo=1,
    #                        log_level=0)
    #     m = nlp1(solver=test_solver)
    #     status = solve(m)
    #
    #     @test status == :Optimal
    #     @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    #     @test m.internalModel.logs[:n_iter] == 1
    # end

    @testset " Validation Test || PBT-AMP-CONV || basic solve || examples/nlp1.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                           bilinear_convexhull=true,
                           monomial_convexhull=true,
                           presolve_bound_tightening=true,
                           presolve_bound_tightening_algo=2,
                           log_level=0)
        m = nlp1(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 1
    end

    @testset " Validation Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           bilinear_convexhull=true,
                           monomial_convexhull=true,
                           presolve_bound_tightening=false,
                           log_level=1)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 9
    end

    # Diabled due to performance issue: need to add limit
    # @testset " Validation Test || BT-AMP-CONV || basic solve || examples/nlp3.jl" begin
    #     test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    #                        mip_solver=CbcSolver(logLevel=0),
    #                        bilinear_convexhull=true,
    #                        monomial_convexhull=true,
    #                        presolve_bound_tightening=true,
    #                        presolve_bound_tightening_algo=1,
    #                        log_level=0)
    #     m = nlp3(solver=test_solver)
    #     status = solve(m)
    #
    #     @test status == :Optimal
    #     @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    #     @test m.internalModel.logs[:n_iter] == 9
    # end

    # Diabled due to solver scheme mismatch (PajaritoSolver doesn't take pure continuous problem)
    # @testset "Validation Test || PBT-AMP-CONV || basic solve || examples/nlp3.jl" begin
    #     test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    #                        mip_solver=CbcSolver(logLevel=0),
    #                        bilinear_convexhull=true,
    #                        monomial_convexhull=true,
    #                        presolve_bound_tightening=true,
    #                        presolve_bound_tightening_algo=2,
    #                        log_level=0)
    #     m = nlp3(solver=test_solver)
    #     status = solve(m)
    #
    #     @test status == :Optimal
    #     @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    #     @test m.internalModel.logs[:n_iter] == 1
    # end

    @testset " Validation Test || AMP || special problem || ... " begin
        test_solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                               discretization_abs_width_tol=1e-2,
                               discretization_ratio=8,
                               maxiter=6,
                               presolve_bound_tightening = false,
                               presolve_bound_tightening_algo = 1,
                               presolve_bt_output_tol = 1e-1,
                               log_level=0)

        m = circle(solver=test_solver)
        solve(m)

        @test isapprox(m.objVal, 1.4142135534556992; atol=1e-3)
    end
end
