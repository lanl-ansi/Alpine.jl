"""
Helper functions
"""

convex(label::Symbol)::Bool = (label == :convex)
concave(label::Symbol)::Bool = (label == :concave)
linear(label::Symbol)::Bool = (label == :linear)

convex(vertex::DAGVertex)::Bool = convex(vertex.convexity)
concave(vertex::DAGVertex)::Bool = concave(vertex.convexity)
linear(vertex::DAGVertex)::Bool = linear(vertex.convexity)

convex(expr::AlpineExpr)::Bool = convex(expr.convexity)
concave(expr::AlpineExpr)::Bool = concave(expr.convexity)
linear(expr::AlpineExpr)::Bool = linear(expr.convexity)


