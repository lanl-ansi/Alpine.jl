
@testset "function convexity detection" begin
    optimizer = nlp_optimizer()
    num_variables = 4
    x = MOI.add_variables(optimizer, num_variables)
    
    quadratic_terms = []
    affine_terms = Alpine.SAF[]
    constant = 0.0
    push!(quadratic_terms, Alpine.SQT(1.0, x[1], x[1]))
    push!(quadratic_terms, Alpine.SQT(1.0, x[2], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :convex

    push!(quadratic_terms, Alpine.SQT(-2.0, x[1], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :convex
    
    quadratic_terms = []
    push!(quadratic_terms, Alpine.SQT(-1.0, x[1], x[1]))
    push!(quadratic_terms, Alpine.SQT(-1.0, x[2], x[2]))
    quadratic_function = Alpine.SQF(affine_terms, quadratic_terms, constant)
    @test Alpine.get_convexity(quadratic_function) == :concave


    
end 

"""
var<> dp("dp",0.5,10.);
    var<int> ip("ip",-1,1);
    auto exp = log(dp) + sqrt(ip);
    exp.print_symbolic();
    CHECK(exp.is_concave());
    var<> p("p");
    var<> q("q");
    auto cc = p*p + q*q;
    cc.print_symbolic();
    CHECK(cc.is_convex());
    auto cc1 = cc * -1;
    cc1.print_symbolic();
    CHECK(cc1.is_concave());
    cc1 += 2*p*q;
    cc1.print_symbolic();
    CHECK(cc1.is_rotated_soc());
    param<int> aa("aa");
    aa = -1;
    aa = -3;
    auto ff = (aa)*p*p;
    ff.print_symbolic();
    CHECK(ff.is_concave());
    ff *= aa;
    ff.print_symbolic();
    CHECK(ff.is_convex());
    ff *= -1;
    ff.print_symbolic();
    CHECK(ff.is_concave());
    ff *= aa;
    ff.print_symbolic();
    CHECK(ff.is_convex());
    ff += aa*(ip + dp)*q*q;
    ff.print_symbolic();
    CHECK(!ff.is_convex());
    CHECK(!ff.is_concave());
    CHECK(ff.is_polynomial());
    CHECK(ff.get_nb_vars()==4);
    param<> b("b");
    b = 1;
    auto fn = p*p + q*q;
    fn.print_symbolic();
    CHECK(fn.is_convex());
    fn -= 2*p*q;
    fn.print_symbolic();
    fn += expo(p);
    fn.print_symbolic();
    CHECK(fn.is_convex());
"""