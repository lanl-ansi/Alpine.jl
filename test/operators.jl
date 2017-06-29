@testset "Expression Operator Tests" begin
    @testset " Validation Test on expression parsing : part1 " begin
        m = Model(solver=PODSolver(nlp_local_solver=IpoptSolver(),
								   mip_solver=CbcSolver(OutputFlag=0),
								   presolve_perform_bound_tightening=true,
								   presolve_bound_tightening_algo=2,
								   presolve_bt_output_tolerance=1e-1,
								   log_level=1))
        @variable(m, x[1:4]>=0)
        @NLconstraint(m, x[1]^2 >= 1)  					# Basic monomial x[5]=x[1]^2

        @NLconstraint(m, x[1]*x[2] <= 1)				# x[6] <= 1 : x[6] = x[1]*x[2]
        @NLconstraint(m, x[1]^2 * x[2]^2 <= 1)          # x[5] + x[7] <= 1 : x[7] = x[2]^2

        @NLconstraint(m, x[1]*(x[2]*x[3]) >= 1)         # x[9] >= 1 : x[8] = x[2] * x[3] && x[9] = x[1]*x[8]
        @NLconstraint(m, x[1]^2*(x[2]^2 * x[3]^3) <= 1) # x[12] <= 1 : x[10] = x[3] ^ 2 && x[11] = x[7] * x[10] && x[12] = x[11]*[5]

        @test 1 == 1

    end

    @testset " Validation Test on expression parsing : part2" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, (x[1]*x[2]) * x[3] >= 1)       # x[7] >= 1 : x[5] = x[1]*x[2], x[6] = x[5]*x[3]
        @NLconstraint(m, (x[1]^2 * x[2]^2) * x[3]^2 <= 1)

        @NLconstraint(m, x[1] * (x[2]^2 * x[3]^2) >= 1)
        @NLconstraint(m, (x[1]^2 * x[2]) * x[3]^2 <= 1)
        @NLconstraint(m, x[1]^2 * (x[2] * x[3]) >= 1)
        @NLconstraint(m, (x[1] * x[2]^2) * x[3] <= 1)

        @test 1 == 1
    end

    @testset "Validation Test on expression parsing: part3" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, ((x[1]*x[2])*x[3])*x[4] >= 1)
        @NLconstraint(m, ((x[1]^2*x[2])*x[3])*x[4] <= 1)
        @NLconstraint(m, ((x[1]*x[2]^2)*x[3])*x[4] >= 1)
        @NLconstraint(m, ((x[1]*x[2])*x[3]^2)*x[4] <= 1)
        @NLconstraint(m, ((x[1]*x[2])*x[3])*x[4]^2 >= 1)
        @NLconstraint(m, ((x[1]^2*x[2]^2)*x[3]^2)*x[4]^2 <=1)

    end

    @testset "Validation Test on expression parsing: part4" begin

        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*(x[2]*(x[3]*x[4])) >= 1)
        @NLconstraint(m, x[1]^2*(x[2]*(x[3]*x[4])) <= 1)
        @NLconstraint(m, x[1]*(x[2]^2*(x[3]*x[4])) >= 1)
        @NLconstraint(m, x[1]*(x[2]*(x[3]^2*x[4])) <= 1)
        @NLconstraint(m, x[1]*(x[2]*(x[3]*x[4]^2)) >= 1)
        @NLconstraint(m, x[1]^2*(x[2]^2*(x[3]^2*x[4]^2)) <= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part5" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*x[2]*x[3] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3] >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3] >= 1)
        @NLconstraint(m, x[1]*x[2]*x[3]^2 >= 1)
        @NLconstraint(m, x[1]^2*x[2]^2*x[3] >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]^2 >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]^2 >= 1)
        @NLconstraint(m, x[1]^2*x[2]^2*x[3]^2 >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part6" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, (x[1]*x[2])*(x[3]*x[4]) >= 1)
        @NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]) >= 1)
        @NLconstraint(m, (x[1]*x[2]^2)*(x[3]*x[4]) >= 1)
        @NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]) >= 1)
        @NLconstraint(m, (x[1]*x[2])*(x[3]^2*x[4]^2) >= 1)
        @NLconstraint(m, (x[1]^2*x[2])*(x[3]*x[4]^2) >= 1)
        @NLconstraint(m, (x[1]^2*x[2])*(x[3]^2*x[4]) >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part7" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*x[2]*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]*x[4] >= 1)
        @NLconstraint(m, x[1]*x[2]*x[3]^2*x[4]^2 >= 1)
        @NLconstraint(m, x[1]*x[2]^2*x[3]^2*x[4] >= 1)
        @NLconstraint(m, x[1]^2*x[2]*x[3]*x[4]^2 >= 1)
        @NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*x[4]^2 >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part8" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, (x[1]*x[2]*x[3])*x[4] >= 1)
        @NLconstraint(m, (x[1]^2*x[2]*x[3])*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2]^2*x[3])*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2]*x[3]^2)*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2]*x[3])*x[4]^2 >= 1)
        @NLconstraint(m, (x[1]^2*x[2]^2*x[3]^2)*x[4]^2 >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part9" begin

        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*(x[2]*x[3])*x[4] >= 1)
        @NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4] >= 1)
        @NLconstraint(m, x[1]*(x[2]^2*x[3])*x[4] >= 1)
        @NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4] >= 1)
        @NLconstraint(m, x[1]*(x[2]*x[3])*x[4]^2 >= 1)
        @NLconstraint(m, x[1]^2*(x[2]*x[3])*x[4]^2 >= 1)
        @NLconstraint(m, x[1]*(x[2]^2*x[3]^2)*x[4] >= 1)
        @NLconstraint(m, x[1]^2*(x[2]^2*x[3])*x[4] >= 1)
        @NLconstraint(m, x[1]*(x[2]*x[3]^2)*x[4]^2 >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part10" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*(x[2]*x[3]*x[4]) >= 1)
        @NLconstraint(m, x[1]^2*(x[2]*x[3]*x[4]) >= 1)
        @NLconstraint(m, x[1]*(x[2]^2*x[3]*x[4]) >= 1)
        @NLconstraint(m, x[1]*(x[2]*x[3]^2*x[4]) >= 1)
        @NLconstraint(m, x[1]*(x[2]*x[3]*x[4]^2) >= 1)
        @NLconstraint(m, x[1]^2*(x[2]*x[3]^2*x[4]) >= 1)
        @NLconstraint(m, x[1]^2*(x[2]^2*x[3]^2*x[4]^2) >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part11" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, x[1]*x[2]*(x[3]*x[4]) >= 1)
        @NLconstraint(m, x[1]^2*x[2]*(x[3]*x[4]) >= 1)
        @NLconstraint(m, x[1]*x[2]^2*(x[3]*x[4]) >= 1)
        @NLconstraint(m, x[1]*x[2]*(x[3]^2*x[4]) >= 1)
        @NLconstraint(m, x[1]*x[2]*(x[3]*x[4]^2) >= 1)
        @NLconstraint(m, x[1]^2*x[2]^2*(x[3]^2*x[4]^2) >= 1)

        @test 1 == 1

    end

    @testset "Validation Test on expression parsing: part12" begin
        m = Model()
        @variable(m, x[1:4]>=0)

        @NLconstraint(m, (x[1]*x[2])*x[3]*x[4] >= 1)
        @NLconstraint(m, (x[1]^2*x[2])*x[3]*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2]^2)*x[3]*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4] >= 1)
        @NLconstraint(m, (x[1]*x[2])*x[3]*x[4]^2 >= 1)
        @NLconstraint(m, (x[1]*x[2])*x[3]^2*x[4]^2 >= 1)
        @NLconstraint(m, (x[1]^2*x[2]^2)*x[3]^2*x[4]^2 >= 1)

    end

end
