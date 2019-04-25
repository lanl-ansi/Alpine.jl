"""
Helper functions
"""

convex(label::Symbol)::Bool = (label == :convex || label == :linear)
concave(label::Symbol)::Bool = (label == :concave || label == :linear)
linear(label::Symbol)::Bool = (label == :linear || (convex(label) && concave(label)))

convex(vertex::DAGVertex)::Bool = convex(vertex.convexity) || linear(vertex.convexity)
concave(vertex::DAGVertex)::Bool = concave(vertex.convexity) || linear(vertex.convexity)
linear(vertex::DAGVertex)::Bool = linear(vertex.convexity) || (convex(vertex) && concave(vertex))

convex(expr::AlpineExpr)::Bool = convex(expr.convexity) || linear(expr.convexity)
concave(expr::AlpineExpr)::Bool = concave(expr.convexity) || linear(expr.convexity)
linear(expr::AlpineExpr)::Bool = linear(expr.convexity) || (convex(expr) && concave(expr))

function get_convexity(is_convex::Bool, is_concave::Bool)::Symbol 
    (is_convex && is_concave) && (return :linear)
    (is_convex) && (return :convex)
    (is_concave) && (return :concave)
    return :undet
end 

negative(interval::Interval{Float64})::Bool = interval.hi < 0
non_negative(interval::Interval{Float64})::Bool = interval.lo >= 0
positive(interval::Interval{Float64})::Bool = interval.lo > 0 
non_positive(interval::Interval{Float64})::Bool = interval.hi <= 0