@testset " Validation Test || AMP-TMC || basic solve || exampls/nlp1.jl" begin

    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=false,
                       monomial_convexhull=false,
                       presolve_bt=false,
                       presolve_bp=true,
                       presolve_bt_output_tol=1e-1,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset " Validation Test || AMP-TMC || basic solve || examples/nlp3.jl (3 iterations)" begin

    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
    					   mip_solver=CbcSolver(logLevel=0),
                           bilinear_convexhull=false,
                           monomial_convexhull=false,
                           presolve_bp=true,
    					   loglevel=100,
                           maxiter=3,
    					   presolve_bt_width_tol=1e-3,
    					   presolve_bt=false,
    					   disc_var_pick=0)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :UserLimits

    @test isapprox(m.objVal, 7049.247897696512; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 3
end

@testset " Validation Test || AMP-TMC || minimum-vertex solving || examples/nlp3.jl (3 iterations)" begin

    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0, max_iter=9999),
                           mip_solver=CbcSolver(logLevel=0),
                           bilinear_convexhull=false,
                           monomial_convexhull=false,
                           presolve_bp=true,
                           disc_var_pick=1,
                           loglevel=100,
                           maxiter=3,
                           presolve_bt_width_tol=1e-3,
                           presolve_bt=false)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :UserLimits
    @test isapprox(m.objVal, 7049.247897696512; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 3
    @test isapprox(m.objBound, 3647.178; atol=1e-2)
end

@testset " Validation Test || PBT-AMP-TMC || basic solve || exampls/nlp1.jl" begin

    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
							   mip_solver=pavito_solver,
                               bilinear_convexhull=false,
                               monomial_convexhull=false,
							   presolve_bt=true,
							   presolve_bt_algo=2,
                               presolve_bp=true,
                               presolve_bt_output_tol=1e-1,
							   loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 2
end

@testset " Validation Test || BT-AMP-TMC || basic solve || examples/nlp3.jl" begin

    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
							   mip_solver=CbcSolver(logLevel=0),
                               bilinear_convexhull=false,
							   loglevel=100,
                               maxiter=3,
							   presolve_bt_width_tol=1e-3,
							   presolve_bt_output_tol=1e-1,
							   presolve_bt=true,
                               presolve_bt_algo=1,
                               presolve_bp=true,
							   presolve_maxiter=2,
                               presolve_track_time=true,
							   disc_var_pick=max_cover_var_picker)
    m = nlp3(solver=test_solver)

    status = solve(m)

    @test status == :UserLimits
    @test m.internalModel.logs[:n_iter] == 3
    @test m.internalModel.logs[:bt_iter] == 2
end

@testset " Validation Test || PBT-AMP-TMC || basic solve || examples/nlp3.jl" begin

    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                               mip_solver=pavito_solver,
                               bilinear_convexhull=false,
                               loglevel=100,
                               maxiter=2,
                               presolve_bt=true,
                               presolve_bt_width_tol=1e-3,
                               presolve_bt_output_tol=1e-1,
                               presolve_bt_algo=2,
                               presolve_bp=true,
                               presolve_maxiter=2,
                               disc_var_pick=max_cover_var_picker)

    m = nlp3(solver=test_solver)
    status = solve(m)
    @test status == :UserLimits
    @test m.internalModel.logs[:n_iter] == 2
    @test m.internalModel.logs[:bt_iter] == 2
end

@testset " Validation Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=true,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset " Validation Test || PBT-AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=true,
                       presolve_bp=true,
                       presolve_bt_algo=2,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 1
end

@testset " Validation Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=CbcSolver(logLevel=0),
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=false,
                       loglevel=100)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 9
end

@testset " Validation Test || AMP || basic solve || examples/circle.jl" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=pavito_solver,
                           disc_abs_width_tol=1e-2,
                           disc_ratio=8,
                           maxiter=6,
                           presolve_bt = false,
                           presolve_bt_algo = 1,
                           presolve_bt_output_tol = 1e-1,
                           loglevel=100)

    m = circle(solver=test_solver)
    solve(m)

    @test isapprox(m.objVal, 1.4142135534556992; atol=1e-3)
end

@testset " Validation Test || AMP || basic solve || examples/circleN.jl" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=pavito_solver,
                           disc_abs_width_tol=1e-2,
                           disc_ratio=8,
                           presolve_bt = false,
                           presolve_bt_algo = 1,
                           presolve_bt_output_tol = 1e-1,
                           loglevel=100)

    m = circleN(solver=test_solver, N=4)
    solve(m)
    @test isapprox(m.objVal, 2.0; atol=1e-3)
end

@testset " Validation Test || AMP-CONV-FACET || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=true,
                       convhull_formulation="facet",
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset " Validation Test || AMP-CONV-MINIB || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=true,
                       convhull_formulation="mini",
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset " Validation Test || AMP || multi4N || N = 2 || exprmode=1:11" begin

    objBoundVec = Any[4.68059, 12.0917, 8.94604, 10.0278, 8.100, 6.6384, 12.5674, 7.3975, 6.0292, 7.9146, 7.8830]
    objValVec = Any[2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    for i in 1:11
        test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(logLevel=0),
                               disc_abs_width_tol=1e-2,
                               maxiter=4,
                               presolve_bp=false,
                               presolve_bt=false,
                               loglevel=1)

        m = multi4N(solver=test_solver, N=2, exprmode=i)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(getobjectivevalue(m), objValVec[i];atol=1e-3)
        @test isapprox(getobjectivebound(m), objBoundVec[i];atol=1e-3)
    end
end

@testset " Validation Test || AMP || multi2 || exprmode=1:11" begin

    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           maxiter=4,
                           presolve_bp=false,
                           presolve_bt=false,
                           loglevel=1)

    m = multi2(solver=test_solver)
    status = solve(m)

    @test status == :UserLimits
    @test isapprox(getobjectivevalue(m), 1.00000;atol=1e-3)
    @test isapprox(getobjectivebound(m), 1.0074;atol=1e-3)
end

@testset " Validation Test || AMP || multi3N || N = 2 || exprmode=1:11" begin

    objBoundVec = Any[2.97186, 3.85492, 4.23375]
    objValVec = Any[2.0, 2.0, 2.0]
    for i in 1:3
        test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                               mip_solver=CbcSolver(logLevel=0),
                               disc_abs_width_tol=1e-2,
                               maxiter=4,
                               presolve_bp=false,
                               presolve_bt=false,
                               loglevel=1)

        m = multi3N(solver=test_solver, N=2, exprmode=i)
        status = solve(m)

        @test status == :UserLimits
        @test isapprox(getobjectivevalue(m), objValVec[i];atol=1e-3)
        @test isapprox(getobjectivebound(m), objBoundVec[i];atol=1e-3)
    end
end

@testset " Validation Test || AMP || multiKND || K = 3, N = 3, D = 0 " begin

    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           maxiter=3,
                           presolve_bp=false,
                           presolve_bt=false,
                           loglevel=1)

    m = multiKND(solver=test_solver, randomub=50, K=3, N=3, D=0)
    status = solve(m)

    @test status == :UserLimits
    @test isapprox(getobjectivevalue(m),3.0000000824779454;atol=1e-3)
    @test isapprox(getobjectivebound(m),12.054604248046875;atol=1e-3)
end

@testset " Validation Test || AMP-CONV-FACET || basic solve || examples/nlp3.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=GLPKSolverMIP(msg_lev=0),
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=false,
                       maxiter=4,
                       convhull_formulation="facet",
                       loglevel=100)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :UserLimits
    @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    @test m.objBound >= 6717.00
    @test m.objBound <= 6718.00
    @test m.internalModel.logs[:n_iter] == 4
end

@testset " Validation Test || AMP-CONV-MINIB || basic solve || examples/nlp3.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=GLPKSolverMIP(msg_lev=0),
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=false,
                       maxiter=4,
                       convhull_formulation="mini",
                       loglevel=100)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :UserLimits
    @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    @test m.objBound >= 6717.00
    @test m.objBound <= 6718.00
    @test m.internalModel.logs[:n_iter] == 4
end

@testset " Validation Test || AMP || DISC-RATIO || examples/nlp3.jl " begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=false,
                           disc_ratio=18,
                           maxiter=1,
                           presolve_bp=true,
                           presolve_bt=false,
                           loglevel=100)

    m = nlp3(solver=test_solver)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 18
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/nlp3.jl " begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=true,
                           presolve_bt=false,
                           loglevel=100)

    m = nlp3(solver=test_solver)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 14
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/castro2m2.jl " begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=true,
                           presolve_bt=false,
                           loglevel=100)

    m = castro2m2(solver=test_solver)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 8
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/multi3N.jl exprmode=2" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=true,
                           presolve_bt=false,
                           loglevel=100)

    m = multi3N(solver=test_solver, N=3, exprmode=1)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 16
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/multi3N.jl exprmode=2" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=false,
                           presolve_bt=false,
                           loglevel=100)

    m = multi3N(solver=test_solver, N=3, exprmode=1)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 20
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/multi4N.jl exprmode=1" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=true,
                           presolve_bt=false,
                           loglevel=100)

    m = multi4N(solver=test_solver, N=2, exprmode=1)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 12
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/multi4N.jl exprmode=2" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=false,
                           presolve_bt=false,
                           loglevel=100)

    m = multi4N(solver=test_solver, N=2, exprmode=1)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 20
end

@testset " Validation Test || AMP || DISC-RATIO-BRANCH || examples/multi4N.jl exprmode=2" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_abs_width_tol=1e-2,
                           disc_ratio_branch=true,
                           maxiter=1,
                           presolve_bp=true,
                           presolve_bt=false,
                           loglevel=100)

    m = multi4N(solver=test_solver, N=2, exprmode=2)
    solve(m)

    @test m.internalModel.logs[:n_iter] == 1
    @test m.internalModel.disc_ratio == 20
end

@testset "Operator :: bmpl && binlin && binprod solve test I" begin
    test_solver=PODSolver(minlp_solver=pavito_solver,
                          nlp_solver=IpoptSolver(print_level=0),
                          mip_solver=CbcSolver(logLevel=0),
                          presolve_bt=false,
                          loglevel=100)

    m = bpml_lnl(test_solver)
    solve(m)
    @test isapprox(m.objVal, 0.3; atol=1e-6)
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[1]), :(x[6])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[2]), :(x[7])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[3]), :(x[8])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[4]), :(x[9])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[5]), :(x[10])])

    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[9])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[10])]][:nonlinear_type] == :BINLIN
end

@testset "Operator :: bmpl && binlin && binprod solve test II" begin
    test_solver=PODSolver(minlp_solver=pavito_solver,
                          nlp_solver=IpoptSolver(print_level=0),
                          mip_solver=CbcSolver(logLevel=0),
                          presolve_bt=false,
                          loglevel=100)

    m = bpml_binl(test_solver)
    solve(m)
    @test isapprox(m.objVal, 15422.058099086951; atol=1e-4)

    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[6]), :(x[7])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[7]), :(x[8])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[8]), :(x[9])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[9]), :(x[10])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[10]), :(x[6])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[1]), :(x[11])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[2]), :(x[13])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[3]), :(x[15])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[4]), :(x[17])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[5]), :(x[19])])

    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[13])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[15])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[17])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[19])]][:nonlinear_type] == :BINLIN
end

@testset "Operator :: bmpl && binlin && binprod solve test II" begin
    test_solver=PODSolver(minlp_solver=pavito_solver,
                          nlp_solver=IpoptSolver(print_level=0),
                          mip_solver=pavito_solver,
                          loglevel=100)

    m = bpml_monl(test_solver)
    solve(m)

    @test isapprox(m.objVal, 19812.920945096557;atol=1e-4)
    @test isapprox(m.objBound, 19812.9205;atol=1e-3)

    nlk1 = Expr[:(x[9]), :(x[9])]
    nlk2 = Expr[:(x[10]), :(x[10])]
    nlk3 = Expr[:(x[8]), :(x[8])]
    nlk4 = Expr[:(x[3]), :(x[15])]
    nlk5 = Expr[:(x[6]), :(x[6])]
    nlk6 = Expr[:(x[2]), :(x[13])]
    nlk7 = Expr[:(x[4]), :(x[17])]
    nlk8 = Expr[:(x[5]), :(x[19])]
    nlk9 = Expr[:(x[7]), :(x[7])]
    nlk10 = Expr[:(x[1]), :(x[11])]

    @test m.internalModel.nonconvex_terms[nlk1][:id] == 7
    @test m.internalModel.nonconvex_terms[nlk2][:id] == 9
    @test m.internalModel.nonconvex_terms[nlk3][:id] == 5
    @test m.internalModel.nonconvex_terms[nlk4][:id] == 6
    @test m.internalModel.nonconvex_terms[nlk5][:id] == 1
    @test m.internalModel.nonconvex_terms[nlk6][:id] == 4
    @test m.internalModel.nonconvex_terms[nlk7][:id] == 8
    @test m.internalModel.nonconvex_terms[nlk8][:id] == 10
    @test m.internalModel.nonconvex_terms[nlk9][:id] == 3
    @test m.internalModel.nonconvex_terms[nlk10][:id] == 2

    @test m.internalModel.nonconvex_terms[nlk1][:y_idx] == 17
    @test m.internalModel.nonconvex_terms[nlk2][:y_idx] == 19
    @test m.internalModel.nonconvex_terms[nlk3][:y_idx] == 15
    @test m.internalModel.nonconvex_terms[nlk4][:y_idx] == 16
    @test m.internalModel.nonconvex_terms[nlk5][:y_idx] == 11
    @test m.internalModel.nonconvex_terms[nlk6][:y_idx] == 14
    @test m.internalModel.nonconvex_terms[nlk7][:y_idx] == 18
    @test m.internalModel.nonconvex_terms[nlk8][:y_idx] == 20
    @test m.internalModel.nonconvex_terms[nlk9][:y_idx] == 13
    @test m.internalModel.nonconvex_terms[nlk10][:y_idx] == 12

    @test m.internalModel.nonconvex_terms[nlk1][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[nlk2][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[nlk3][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[nlk4][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[nlk5][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[nlk6][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[nlk7][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[nlk8][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[nlk9][:nonlinear_type] == :MONOMIAL
    @test m.internalModel.nonconvex_terms[nlk10][:nonlinear_type] == :BINLIN
end

@testset "Operator :: bmpl && binlin && binprod solve test III with negative bounds" begin
    test_solver=PODSolver(minlp_solver=pavito_solver,
                          nlp_solver=IpoptSolver(print_level=0),
                          mip_solver=pavito_solver,
                          disc_var_pick=1,
                          loglevel=100)

    m = bpml_negative(test_solver)
    solve(m)

    @test m.objVal <= 12789.0
    @test m.objBound >= 12788.0

    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[6]), :(x[7])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[7]), :(x[8])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[8]), :(x[9])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[9]), :(x[10])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[10]), :(x[6])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[1]), :(x[11])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[2]), :(x[13])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[3]), :(x[15])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[4]), :(x[17])])
    @test haskey(m.internalModel.nonconvex_terms, Expr[:(x[5]), :(x[19])])

    @test m.internalModel.nonconvex_terms[Expr[:(x[6]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[7]), :(x[8])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[8]), :(x[9])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[11])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[13])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[15])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[4]), :(x[17])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[5]), :(x[19])]][:nonlinear_type] == :BINLIN
end

@testset "Embedding Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=true,
                       convhull_ebd=true,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset "Embedding Test || PBT-AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=true,
                       presolve_bp=true,
                       presolve_bt_algo=2,
                       convhull_ebd=true,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 1
end

@testset "Embedding Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=CbcSolver(logLevel=0),
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=false,
                       convhull_ebd=true,
                       loglevel=100)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 9
end

@testset "Embedding Test || AMP || special problem || ... " begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=pavito_solver,
                           disc_abs_width_tol=1e-2,
                           disc_ratio=8,
                           maxiter=6,
                           presolve_bt=false,
                           presolve_bp=true,
                           presolve_bt_algo=1,
                           presolve_bt_output_tol=1e-1,
                           convhull_ebd=true,
                           loglevel=100)

    m = circle(solver=test_solver)
    solve(m)
    @test isapprox(m.objVal, 1.4142135534556992; atol=1e-3)
end

@testset "Embedding IBS Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=true,
                       convhull_ebd=true,
                       convhull_ebd_ibs=true,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset "Embedding IBS Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=CbcSolver(logLevel=0),
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=false,
                       convhull_ebd=true,
                       convhull_ebd_ibs=true,
                       loglevel=100)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 9
end

@testset "Embedding IBS Test || AMP || special problem || ... " begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=pavito_solver,
                           disc_abs_width_tol=1e-2,
                           disc_ratio=8,
                           maxiter=6,
                           presolve_bt=false,
                           presolve_bp=true,
                           presolve_bt_algo=1,
                           presolve_bt_output_tol=1e-1,
                           convhull_ebd=true,
                           convhull_ebd_ibs=true,
                           loglevel=100)

    m = circle(solver=test_solver)
    solve(m)
    @test isapprox(m.objVal, 1.4142135534556992; atol=1e-3)
end

@testset "Embedding LINK Test || AMP-CONV || basic solve || examples/nlp1.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=pavito_solver,
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=true,
                       convhull_ebd=true,
                       convhull_ebd_link=true,
                       loglevel=100)
    m = nlp1(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 58.38367169858795; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 7
end

@testset "Embedding LINK Test || AMP-CONV || basic solve || examples/nlp3.jl" begin
    test_solver = PODSolver(nlp_solver=IpoptSolver(print_level=0),
                       mip_solver=CbcSolver(logLevel=0),
                       bilinear_convexhull=true,
                       monomial_convexhull=true,
                       presolve_bt=false,
                       presolve_bp=false,
                       convhull_ebd=true,
                       convhull_ebd_link=true,
                       loglevel=100)
    m = nlp3(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 7049.247897696188; atol=1e-4)
    @test m.internalModel.logs[:n_iter] == 9
end

@testset "Embedding LINK Test || AMP || special problem || ... " begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=pavito_solver,
                           disc_abs_width_tol=1e-2,
                           disc_ratio=8,
                           maxiter=6,
                           presolve_bt=false,
                           presolve_bp=true,
                           presolve_bt_algo=1,
                           presolve_bt_output_tol=1e-1,
                           convhull_ebd=true,
                           convhull_ebd_link=true,
                           loglevel=100)

    m = circle(solver=test_solver)
    solve(m)
    @test isapprox(m.objVal, 1.4142135534556992; atol=1e-3)
end

@testset "Algorithm Logic Test || castro4m2 || 1 iteration || Error case" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           maxiter=1,
                           colorful_pod="warmer",
                           presolve_bt=false,
                           loglevel=100)

    m = castro4m2(solver=test_solver)
    status = solve(m)
    @test status == :UserLimits
    @test m.internalModel.status[:local_solve] == :Error
end

@testset " Algorithm Logic Test || blend029_gl || 3 iterations || Infeasible Case" begin

    test_solver=PODSolver(minlp_solver=pavito_solver,
                          nlp_solver=IpoptSolver(print_level=0),
                          mip_solver=CbcSolver(logLevel=0),
                          presolve_bp=true,
                          disc_var_pick=1,
                          loglevel=100,
                          maxiter=3,
                          presolve_bt_width_tol=1e-3,
                          presolve_bt=false)
    m = blend029_gl(solver=test_solver)
    status = solve(m)

    @test status == :UserLimits
    @test m.internalModel.logs[:n_iter] == 3
    @test getobjbound(m) <= 14.0074
end

@testset "Convex Model Solve" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=pavito_solver,
                           maxiter=1,
                           colorful_pod="solarized",
                           presolve_bt=false,
                           loglevel=100)
    m = convex_solve(solver=test_solver)
    status = solve(m)
    @test status == :Optimal
end

@testset "Uniform partitioning" begin
    test_solver=PODSolver(nlp_solver=IpoptSolver(print_level=0),
                           mip_solver=CbcSolver(logLevel=0),
                           disc_add_partition_method = "uniform",
                           disc_uniform_rate = 10,
                           maxiter=1,
                           presolve_bt=false,
                           colorful_pod="random",
                           timeout=100000,
                           loglevel=100)
    m = nlp3(solver=test_solver)
    solve(m)
    @test isapprox(m.objBound, 6561.7156;atol=1e-3)
end

@testset "Algorithm Test with binprod terms" begin

    test_solver = PODSolver(minlp_solver=pavito_solver,
                            nlp_solver=IpoptSolver(print_level=0),
                            mip_solver=CbcSolver(logLevel=0),
                            bilinear_convexhull=true,
                            monomial_convexhull=true,
                            presolve_bp=true,
                            presolve_bt=false,
                            loglevel=100)
    m = binprod_nlp3(solver=test_solver)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(m.objVal, 3651.020370626844;atol=1e-3)
    @test isapprox(m.objBound, 3650.791316892635;atol=1e-3)

    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:y_idx] == 19
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:id] == 6
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:lifted_constr_ref] == :(x[19] == x[2] * x[4])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[4])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:y_idx] == 25
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:id] == 12
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:lifted_constr_ref] == :(x[25] == x[3] * x[8])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[8])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:y_idx] == 22
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:id] == 9
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:lifted_constr_ref] == :(x[22] == x[3] * x[5])
    @test m.internalModel.nonconvex_terms[Expr[:(x[3]), :(x[5])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:y_idx] == 20
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:id] == 7
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:lifted_constr_ref] == :(x[20] == x[12] * x[19])
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[19])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:y_idx] == 14
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:id] == 1
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:lifted_constr_ref] == :(x[14] == x[9] * x[4])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[4])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:y_idx] == 17
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:id] == 4
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:lifted_constr_ref] == :(x[17] == x[1] * x[6])
    @test m.internalModel.nonconvex_terms[Expr[:(x[1]), :(x[6])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:y_idx] == 16
    @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:id] == 3
    @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:lifted_constr_ref] == :(x[16] == x[11] * x[5])
    @test m.internalModel.nonconvex_terms[Expr[:(x[11]), :(x[5])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:y_idx] == 30
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:id] == 17
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:lifted_constr_ref] == :(x[30] == x[10] * x[11])
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[11])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:y_idx] == 29
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:id] == 16
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:lifted_constr_ref] == :(x[29] == x[9] * x[10] * x[11])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[10]), :(x[11])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:y_idx] == 21
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:id] == 8
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:lifted_constr_ref] == :(x[21] == x[2] * x[7])
    @test m.internalModel.nonconvex_terms[Expr[:(x[2]), :(x[7])]][:nonlinear_type] == :BILINEAR
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:y_idx] == 32
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:id] == 19
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:lifted_constr_ref] == :(x[32] == x[9] * x[13])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[13])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:y_idx] == 15
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:id] == 2
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:lifted_constr_ref] == :(x[15] == x[10] * x[6])
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[6])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:y_idx] == 33
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:id] == 20
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:lifted_constr_ref] == :(x[33] == x[10] * x[12])
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[12])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:y_idx] == 28
    @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:id] == 15
    @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:lifted_constr_ref] == :(x[28] == x[27] * x[5])
    @test m.internalModel.nonconvex_terms[Expr[:(x[27]), :(x[5])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:y_idx] == 26
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:id] == 13
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:lifted_constr_ref] == :(x[26] == x[13] * x[25])
    @test m.internalModel.nonconvex_terms[Expr[:(x[13]), :(x[25])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:y_idx] == 27
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:id] == 14
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:lifted_constr_ref] == :(x[27] == x[9] * x[12])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[12])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:y_idx] == 23
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:id] == 10
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:lifted_constr_ref] == :(x[23] == x[10] * x[13])
    @test m.internalModel.nonconvex_terms[Expr[:(x[10]), :(x[13])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:y_idx] == 24
    @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:id] == 11
    @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:lifted_constr_ref] == :(x[24] == x[23] * x[22])
    @test m.internalModel.nonconvex_terms[Expr[:(x[23]), :(x[22])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:y_idx] == 31
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:id] == 18
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:lifted_constr_ref] == :(x[31] == x[12] * x[13])
    @test m.internalModel.nonconvex_terms[Expr[:(x[12]), :(x[13])]][:nonlinear_type] == :BINPROD
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:y_idx] == 18
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:id] == 5
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:lifted_constr_ref] == :(x[18] == x[9] * x[17])
    @test m.internalModel.nonconvex_terms[Expr[:(x[9]), :(x[17])]][:nonlinear_type] == :BINLIN
    @test m.internalModel.bounding_constr_mip[1][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[1][:vars] == Any[:(x[14]), :(x[15])]
    @test m.internalModel.bounding_constr_mip[1][:coefs] == Any[0.0025, 0.0025]
    @test m.internalModel.bounding_constr_mip[1][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[1][:cnt] == 2
    @test m.internalModel.bounding_constr_mip[2][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[2][:vars] == Any[:(x[14]), :(x[5]), :(x[7])]
    @test m.internalModel.bounding_constr_mip[2][:coefs] == Any[-0.0025, 0.0025, 0.0025]
    @test m.internalModel.bounding_constr_mip[2][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[2][:cnt] == 3
    @test m.internalModel.bounding_constr_mip[3][:rhs] == 1.0
    @test m.internalModel.bounding_constr_mip[3][:vars] == Any[:(x[16]), :(x[8])]
    @test m.internalModel.bounding_constr_mip[3][:coefs] == Any[-0.01, 0.01]
    @test m.internalModel.bounding_constr_mip[3][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[3][:cnt] == 2
    @test m.internalModel.bounding_constr_mip[4][:rhs] == 83333.333
    @test m.internalModel.bounding_constr_mip[4][:vars] == Any[:(x[1]), :(x[18]), :(x[14])]
    @test m.internalModel.bounding_constr_mip[4][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[4][:cnt] == 3
    @test m.internalModel.bounding_constr_mip[5][:rhs] == 0.0
    @test m.internalModel.bounding_constr_mip[5][:vars] == Any[:(x[20]), :(x[21]), :(x[4]), :(x[5])]
    @test m.internalModel.bounding_constr_mip[5][:coefs] == Any[1.0, -1.0, -1250.0, 1250.0]
    @test m.internalModel.bounding_constr_mip[5][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[5][:cnt] == 4
    @test m.internalModel.bounding_constr_mip[6][:rhs] == -1.25e6
    @test m.internalModel.bounding_constr_mip[6][:vars] == Any[:(x[24]), :(x[26]), :(x[28])]
    @test m.internalModel.bounding_constr_mip[6][:coefs] == Any[1.0, -1.0, -2500.0]
    @test m.internalModel.bounding_constr_mip[6][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[6][:cnt] == 3
    @test m.internalModel.bounding_constr_mip[7][:rhs] == 0.0
    @test m.internalModel.bounding_constr_mip[7][:vars] == Any[:(x[29])]
    @test m.internalModel.bounding_constr_mip[7][:coefs] == Any[1.0]
    @test m.internalModel.bounding_constr_mip[7][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[7][:cnt] == 1
    @test m.internalModel.bounding_constr_mip[8][:rhs] == 0.0
    @test m.internalModel.bounding_constr_mip[8][:vars] == Any[:(x[30]), :(x[31])]
    @test m.internalModel.bounding_constr_mip[8][:coefs] == Any[1.0, -1.0]
    @test m.internalModel.bounding_constr_mip[8][:sense] == :(>=)
    @test m.internalModel.bounding_constr_mip[8][:cnt] == 2
    @test m.internalModel.bounding_constr_mip[9][:rhs] == 0.0
    @test m.internalModel.bounding_constr_mip[9][:vars] == Any[:(x[32]), :(x[33])]
    @test m.internalModel.bounding_constr_mip[9][:coefs] == Any[1.0, -1.0]
    @test m.internalModel.bounding_constr_mip[9][:sense] == :(<=)
    @test m.internalModel.bounding_constr_mip[9][:cnt] == 2
end
