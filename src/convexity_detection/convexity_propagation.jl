"""

"""
function dag_convexity_propagation!(model::MOI.AbstractOptimizer)

    dag = model.inner.expression_graph 
    dag_lookup = model.inner.dag_lookup
    for d in 1:dag.max_depth 
        for dag_vertex in dag.vertices[d]
            operation = dag_vertex.vertex.args[1]
            children = dag_vertex.children 
            children_interval = Interval{Float64}[]
            children_convexity = Symbol[]
            children_type = Symbol[]
            for child in children
                depth, position = dag_lookup[child]
                push!(children_interval, dag.vertices[depth][position].interval)
                push!(children_convexity, dag.vertices[depth][position].convexity)
                push!(children_type, dag.vertices[depth][position].vertex_type)
            end 
            convexity = get_convexity(operation, 
                children_interval, children_convexity, children_type)
        end 
    end 
    return
end

"""
    function get_convexity(operation::Symbol, 
        children_interval::Vector{Interval{Float64}}, 
        children_convexity::Vector{Symbol},
        children_type::Vector{Symbol})::Symbol
This function dectects convexity of compositions. For details 
see supplementary material of the paper entitled "SUSPECT: 
MINLP Special Structure Detector for Pyomo"
"""
function get_convexity(operation::Symbol, 
    children_interval::Vector{Interval{Float64}}, 
    children_convexity::Vector{Symbol},
    children_type::Vector{Symbol})::Symbol
    
    (:undet in children_convexity) && (return :undet)
    all_constant = all(y -> y == first(children_type), children_type)
    (all_constant) && (return :linear)
    # a = -b
    if operation == :- && length(children_interval) == 1
        is_convex = concave(children_convexity[1])
        is_concave = convex(children_convexity[1])
        return get_convexity(is_convex, is_concave)
    # a = b - c 
    elseif operation == :- && length(children_interval) == 2
        is_convex = convex(children_convexity[1]) && concave(children_convexity[2])
        is_concave = concave(children_convexity[1]) && convex(children_convexity[2])
        return get_convexity(is_convex, is_concave)
    # k-nary addition 
    elseif operation == :+ 
        all_same = all(y -> y == first(children_convexity), children_convexity)
        is_convex = false 
        is_concave = false 
        (!all_same) && (return get_convexity(is_convex, is_concave))
        (all_same && first(children_convexity) == :convex) && (is_convex = true)
        (all_same && first(children_convexity) == :concave) && (is_concave = true)
        return get_convexity(is_convex, is_concave)
    # k-nary multiplication
    elseif operation == :* 
        constant_count = count(y -> y == :constant, children_type)
        if constant_count != length(children_type)-1
            return :undet 
        else 
            child_type = filter(y -> y != :constant, children_type)
            index = findall(y -> y == first(child_type), children_type)
            constant_value = 1.0
            for i in 1:length(children_type)
                (i == index) && (continue)
                constant_value *= children_interval[i].lo
            end
            is_convex = convex(children_convexity[index]) && (constant_value >= 0)
            is_convex |= (concave(children_convexity[index]) && (constant_value <= 0))
            is_concave = concave(children_convexity[index]) && (constant_value >= 0)
            is_concave |= (convex(children_convexity[index]) && (constant_value <= 0))
            return get_convexity(is_convex, is_concave)
        end 
    # a = b/c
    elseif operation == :/
        if children_type[1] != :constant && children_type[2] != :constant 
            return :undet 
        elseif children_type[1] == :constant 
            constant_value = children_interval[1].lo 
            is_convex = (constant_value >= 0) && 
                (concave(children_convexity[2]) && positive(children_interval[2])) 
            is_convex |= (
                (constant_value <= 0) && 
                (convex(children_convexity[2]) && negative(children_interval[2]))
            )
            is_concave = (constant_value >= 0) && 
                (convex(children_convexity[2]) && negative(children_interval[2]))
            is_concave |= (
                (constant_value >= 0) && 
                (concave(children_convexity[2]) && positive(children_interval[2]))
            ) 
            return get_convexity(is_convex, is_concave)
        else
            constant_value = children_interval[2].lo 
            is_convex = convex(children_convexity[1]) && constant_value > 0
            is_convex = is_convex || (
                concave(children_convexity[1]) && constant_value < 0
            )
            is_concave = convex(children_convexity[1]) && constant_value < 0
            is_concave = is_concave || (
                concave(children_convexity[1]) && constant_value > 0
            )
            return get_convexity(is_convex, is_concave)
        end 
    elseif operation == :^ && children_type[1] == :constant 
        constant_value = children_interval[1].lo
        is_convex = (concave(children_convexity[2]) && (0 < constant_value < 1)) || 
            (convex(children_convexity[2]) && constant_value >= 1)
        is_concave = false 
        return get_convexity(is_convex, is_concave)
    elseif operation == :^ && children_type[2] == :constant 
        constant_value = children_interval[2].lo 
        is_convex = false
        is_concave = false
        (constant_value == 1) && (return children_convexity[1])
        is_integer = try isa(Int(constant_value), Int) catch y false end 
        if is_integer 
            if iseven(Int(constant_value)) && constant_value > 0
                (linear(children_convexity[1])) && (return :linear)
                is_convex |= (convex(children_convexity[1]) && non_negative(children_interval[1]))
                is_convex |= (concave(children_convexity[1]) && non_positive(children_interval[1]))
            elseif iseven(Int(constant_value)) && constant_value < 0
                is_convex |= (convex(children_convexity[1]) && non_positive(children_interval[1]))
                is_convex |= (concave(children_convexity[1]) && non_positive(children_interval[1]))
                is_concave |= (convex(children_convexity[1]) && non_negative(children_interval[1]))
                is_concave |= (concave(children_convexity[1]) && non_positive(children_interval[1]))
            elseif isodd(Int(constant_value)) && constant_value > 0
                is_convex |= (convex(children_convexity[1]) && non_negative(children_interval[1]))
                is_concave |= (concave(children_convexity[1]) && non_positive(children_interval[1]))
            end 
        else 
            if constant_value > 1
                is_convex |= (convex(children_convexity[1]) && non_negative(children_interval[1]))
            elseif constant_value < 0 
                is_convex |= (concave(children_convexity[1]) && non_negative(children_interval[1]))
                is_concave |= (convex(children_convexity[1]) && non_negative(children_interval[1]))
            elseif 0 < constant_value < 1
                is_concave |= (concave(children_convexity[1]) && non_negative(children_interval[1]))
            end 
        end
        return get_convexity(is_convex, is_concave)
    elseif operation == :^ && :constant ∉ children_type 
        return :undet 
    elseif operation == :abs
        is_convex = linear(children_convexity[1])
        is_convex |= (convex(children_convexity[1]) && non_negative(children_interval[1]))
        is_convex |= (concave(children_convexity[1]) && non_positive(children_interval[1]))
        is_concave = false
        return get_convexity(is_convex, is_concave)
    elseif operation == :sqrt 
        is_concave = concave(children_convexity[1]) && non_negative(children_interval[1])
        is_convex = false 
        return get_convexity(is_convex, is_concave)
    elseif operation == :exp 
        return get_convexity(convex(children_convexity[1]), false)
    elseif operation == :log
        is_concave = concave(children_convexity[1]) && non_negative(children_interval[1])
        is_convex = false 
        return get_convexity(is_convex, is_concave)
    elseif operation == :sin 
        is_convex = false 
        is_concave = false 
        is_concave |= (linear(children_convexity[1]) && non_negative(sin(children_interval[1])))
        is_concave |= (convex(children_convexity[1]) && 
            non_negative(sin(children_interval[1])) && 
            non_positive(cos(children_interval[1]))
        )
        is_concave |= (concave(children_convexity[1]) && 
            non_negative(sin(children_interval[1])) && 
            non_negative(cos(children_interval[1]))
        )
        is_convex |= (linear(children_convexity[1]) && non_positive(sin(children_interval[1])))
        is_convex |= (concave(children_convexity[1]) && 
            non_positive(sin(children_interval[1])) && 
            non_positive(cos(children_interval[1]))
        )
        is_convex |= (convex(children_convexity[1]) && 
            non_positive(sin(children_interval[1])) && 
            non_negative(cos(children_interval[1]))
        )
        return get_convexity(is_convex, is_concave)
    elseif operation == :cos 
        is_convex = false 
        is_concave = false 
        is_concave |= (linear(children_convexity[1]) && non_negative(cos(children_interval[1])))
        is_concave |= (convex(children_convexity[1]) && 
            non_negative(cos(children_interval[1])) && 
            non_positive(sin(children_interval[1]))
        )
        is_concave |= (concave(children_convexity[1]) && 
            non_negative(cos(children_interval[1])) && 
            non_negative(sin(children_interval[1]))
        )
        is_convex |= (linear(children_convexity[1]) && non_positive(cos(children_interval[1])))
        is_convex |= (concave(children_convexity[1]) && 
            non_positive(cos(children_interval[1])) && 
            non_positive(sin(children_interval[1]))
        )
        is_convex |= (convex(children_convexity[1]) && 
            non_positive(cos(children_interval[1])) && 
            non_negative(sin(children_interval[1]))
        )
        return get_convexity(is_convex, is_concave)
    elseif operation == :tan 
        is_convex = (convex(children_convexity[1]) && non_negative(tan(children_interval[1])))
        is_concave = (concave(children_convexity[1]) && non_positive(tan(children_interval[1])))
    else 
        return :undet 
    end
end 