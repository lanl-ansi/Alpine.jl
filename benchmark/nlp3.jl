benchmark_start = time()
include("../examples/nlp3.jl")
m = nlp3()
status = solve(m)
benchmark_time = time() - benchmark_start
println("Algorithm Final Status: $status")
println("Total time: $(benchmark_time)s")
