#
#   Plan for the following non-linear subtree structure
#   1. bilinear tree
#           (:*)        x*y
#           /  \
#         y     x
#   2. monomial tree
#           (:^)        x^2
#           /  \
#         x     2
#   3. power tree       x^a, a>2
#           (:^)
#           /  \
#         x    >2
#   3. multilinear tree  x*z*y
#           (:*)
#         / | .. \
#       x   z     y
#   4. hierarchical bilinear tree (upon 1)  (x*y)*z
#           (:*)
#           /  \
#         (:*)  z
#         /  \
#        x    y
#   5. sin tree  (can fit 1, 2, 3)   sin(x)
#          (:sin)
#             |
#             x
#   6. cos tree (can fit 1, 2, 3)    cos(x)
#           (:cos)
#             |
#             x
#   7. user-defined tree (can over-ride built-in structural trees)
#
#

"""
    process_expr(expr; kwargs...)
    
High-level warpper for processing expression with sub-tree operators
"""
function process_expr(expr;kwargs...)
    options = Dict(kwargs)

    return
end

function is_bilinear(expr)
    return
end

function is_monomial(expr)
    return
end

function is_multilinear(expr)
    return
end

function is_sin(expr)
    return
end

function is_cos(expr)
    return
end
