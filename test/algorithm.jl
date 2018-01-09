@testset "Solving algorithm tests" begin

    @testset " Validation Test || AMP-TMC || basic solve || exampls/nlp1.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                           bilinear_convexhull=false,
                           monomial_convexhull=false,
                           presolve_bound_tightening=false,
                           bound_basic_propagation=true,
                           presolve_bt_output_tol=1e-1,
                           log_level=100)
        m = nlp1(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 7
    end

    @testset " Validation Test || AMP-TMC || basic solve || examples/nlp3.jl (3 iterations)" begin

        test_solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
        					   mip_solver=CbcSolver(logLevel=0),
                               bilinear_convexhull=false,
                               monomial_convexhull=false,
                               bound_basic_propagation=true,
        					   log_level=1,
                               max_iter=3,
        					   presolve_bt_width_tol=1e-3,
        					   presolve_bound_tightening=false,
        					   discretization_var_pick_algo=0)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :UserLimits

        @test isapprox(m.objVal, 7049.247897696512; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 3
    end

    @testset " Validation Test || PBT-AMP-TMC || basic solve || exampls/nlp1.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    							   mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                                   bilinear_convexhull=false,
                                   monomial_convexhull=false,
    							   presolve_bound_tightening=true,
    							   presolve_bound_tightening_algo=2,
                                   bound_basic_propagation=true,
                                   presolve_bt_output_tol=1e-1,
    							   log_level=0)
        m = nlp1(solver=test_solver)
        status = solve(m)

        # @show m.internalModel.l_var_tight
        # @show m.internalModel.u_var_tight

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 2

    end

    @testset " Validation Test || BT-AMP-TMC || basic solve || examples/nlp3.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
    							   mip_solver=CbcSolver(logLevel=0),
                                   bilinear_convexhull=false,
    							   log_level=0,
                                   max_iter=3,
    							   presolve_bt_width_tol=1e-3,
    							   presolve_bt_output_tol=1e-1,
    							   presolve_bound_tightening=true,
                                   presolve_bound_tightening_algo=1,
                                   bound_basic_propagation=true,
    							   presolve_max_iter=2,
    							   discretization_var_pick_algo=max_cover_var_picker)
        m = nlp3(solver=test_solver)

        status = solve(m)

        @test status == :UserLimits
        @test m.internalModel.logs[:n_iter] == 3
        @test m.internalModel.logs[:bt_iter] == 2

    end

    @testset " Validation Test || PBT-AMP-TMC || basic solve || examples/nlp3.jl" begin

        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                                   mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(), log_level=0),
                                   bilinear_convexhull=false,
                                   log_level=0,
                                   max_iter=2,
                                   presolve_bound_tightening=true,
                                   presolve_bt_width_tol=1e-3,
                                   presolve_bt_output_tol=1e-1,
                                   presolve_bound_tightening_algo=2,
                                   bound_basic_propagation=true,
                                   presolve_max_iter=2,
                                   discretization_var_pick_algo=max_cover_var_picker)

        m = nlp3(solver=test_solver)
        status = solve(m)
        @test status == :UserLimits
        @test m.internalModel.logs[:n_iter] == 2
        @test m.internalModel.logs[:bt_iter] == 2

    end

    @testset " Validation Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(), log_level=0),
                           bilinear_convexhull=true,
                           monomial_convexhull=true,
                           presolve_bound_tightening=false,
                           bound_basic_propagation=true,
                           log_level=0)
        m = nlp1(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 7
    end

    @testset " Validation Test || PBT-AMP-CONV || basic solve || examples/nlp1.jl" begin
        test_solver = PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                           mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                           bilinear_convexhull=true,
                           monomial_convexhull=true,
                           presolve_bound_tightening=true,
                           bound_basic_propagation=true,
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
                           bound_basic_propagation=false,
                           log_level=0)
        m = nlp3(solver=test_solver)
        status = solve(m)

        @test status == :Optimal
        @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
        @test m.internalModel.logs[:n_iter] == 9
    end

    @testset " Validation Test || AMP || special problem || ... " begin
        test_solver=PODSolver(nlp_local_solver=IpoptSolver(print_level=0),
                               mip_solver=PajaritoSolver(cont_solver=IpoptSolver(print_level=0), mip_solver=CbcSolver(logLevel=0), log_level=0),
                               discretization_abs_width_tol=1e-2,
                               discretization_ratio=8,
                               max_iter=6,
                               presolve_bound_tightening=false,
                               bound_basic_propagation=true,
                               presolve_bound_tightening_algo=1,
                               presolve_bt_output_tol=1e-1,
                               log_level=0)

        m = circle(solver=test_solver)
        solve(m)

        @test isapprox(m.objVal, 1.4142135534556992; atol=1e-3)
    end
end

@testset "Solving algorithm tests :: Bin-Lin Solves" begin
    @testset "Operator :: bmpl && linbin && binprod solve test I" begin
        test_solver=PODSolver(minlp_local_solver=PajaritoSolver(mip_solver=CbcSolver(),
                                                              cont_solver=IpoptSolver(),
                                                              log_level=0),
                              nlp_local_solver=IpoptSolver(),
                              mip_solver=CbcSolver(),
                              log_level=100)

        m = bpml_lnl(test_solver)
        solve(m)
        @test isapprox(m.objVal, 0.3; atol=1e-6)
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[1]), :(x[6])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[2]), :(x[7])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[3]), :(x[8])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[4]), :(x[9])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[5]), :(x[10])])

        @test m.internalModel.nonlinear_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[2]), :(x[7])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[3]), :(x[8])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[4]), :(x[9])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[5]), :(x[10])]][:nonlinear_type] == :binlin
    end

    @testset "Operator :: bmpl && linbin && binprod solve test II" begin
        test_solver=PODSolver(minlp_local_solver=PajaritoSolver(mip_solver=CbcSolver(),
                                                              cont_solver=IpoptSolver(),
                                                              log_level=0),
                              nlp_local_solver=IpoptSolver(),
                              mip_solver=CbcSolver(),
                              log_level=100)

        m = bpml_binl(test_solver)
        solve(m)
        @test isapprox(m.objVal, 15422.058099086951; atol=1e-4)

        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[6]), :(x[7])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[7]), :(x[8])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[8]), :(x[9])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[9]), :(x[10])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[10]), :(x[6])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[11]), :(x[1])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[13]), :(x[2])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[15]), :(x[3])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[17]), :(x[4])])
        @test haskey(m.internalModel.nonlinear_terms, Expr[:(x[19]), :(x[5])])

        @test m.internalModel.nonlinear_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[7]), :(x[8])]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[8]), :(x[9])]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[10]), :(x[6])]][:nonlinear_type] == :bilinear
        @test m.internalModel.nonlinear_terms[Expr[:(x[11]), :(x[1])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[13]), :(x[2])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[15]), :(x[3])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[17]), :(x[4])]][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[Expr[:(x[19]), :(x[5])]][:nonlinear_type] == :binlin
    end

    @testset "Operator :: bmpl && linbin && binprod solve test II" begin
        test_solver=PODSolver(minlp_local_solver=PajaritoSolver(mip_solver=CbcSolver(),
                                                              cont_solver=IpoptSolver(),
                                                              log_level=0),
                              nlp_local_solver=IpoptSolver(),
                              mip_solver=PajaritoSolver(mip_solver=CbcSolver(),
                                                        cont_solver=IpoptSolver(),
                                                        log_level=0),
                              discretization_var_pick_algo=1,
                              log_level=100)

        m = bpml_monl(test_solver)
        solve(m)

        @test isapprox(m.objVal, 19812.920945096557;atol=1e-4)
        @test isapprox(m.objBound, 19812.9205;atol=1e-3)

        nlk1 = Expr[:(x[9]), :(x[9])]
        nlk2 = Expr[:(x[10]), :(x[10])]
        nlk3 = Expr[:(x[8]), :(x[8])]
        nlk4 = Expr[:(x[15]), :(x[3])]
        nlk5 = Expr[:(x[6]), :(x[6])]
        nlk6 = Expr[:(x[13]), :(x[2])]
        nlk7 = Expr[:(x[17]), :(x[4])]
        nlk8 = Expr[:(x[19]), :(x[5])]
        nlk9 = Expr[:(x[7]), :(x[7])]
        nlk10 = Expr[:(x[11]), :(x[1])]

        @test m.internalModel.nonlinear_terms[nlk1][:id] == 7
        @test m.internalModel.nonlinear_terms[nlk2][:id] == 9
        @test m.internalModel.nonlinear_terms[nlk3][:id] == 5
        @test m.internalModel.nonlinear_terms[nlk4][:id] == 6
        @test m.internalModel.nonlinear_terms[nlk5][:id] == 1
        @test m.internalModel.nonlinear_terms[nlk6][:id] == 4
        @test m.internalModel.nonlinear_terms[nlk7][:id] == 8
        @test m.internalModel.nonlinear_terms[nlk8][:id] == 10
        @test m.internalModel.nonlinear_terms[nlk9][:id] == 3
        @test m.internalModel.nonlinear_terms[nlk10][:id] == 2

        @test m.internalModel.nonlinear_terms[nlk1][:y_idx] == 17
        @test m.internalModel.nonlinear_terms[nlk2][:y_idx] == 19
        @test m.internalModel.nonlinear_terms[nlk3][:y_idx] == 15
        @test m.internalModel.nonlinear_terms[nlk4][:y_idx] == 16
        @test m.internalModel.nonlinear_terms[nlk5][:y_idx] == 11
        @test m.internalModel.nonlinear_terms[nlk6][:y_idx] == 14
        @test m.internalModel.nonlinear_terms[nlk7][:y_idx] == 18
        @test m.internalModel.nonlinear_terms[nlk8][:y_idx] == 20
        @test m.internalModel.nonlinear_terms[nlk9][:y_idx] == 13
        @test m.internalModel.nonlinear_terms[nlk10][:y_idx] == 12

        @test m.internalModel.nonlinear_terms[nlk1][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk2][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk3][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk4][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk5][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk6][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk7][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk8][:nonlinear_type] == :binlin
        @test m.internalModel.nonlinear_terms[nlk9][:nonlinear_type] == :monomial
        @test m.internalModel.nonlinear_terms[nlk10][:nonlinear_type] == :binlin

    end
end
