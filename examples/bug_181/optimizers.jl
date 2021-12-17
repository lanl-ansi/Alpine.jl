using Juniper, CPLEX, Ipopt, Gurobi

function get_gurobi()
    return optimizer_with_attributes(Gurobi.Optimizer, 
                                    MOI.Silent() => true, 
                                    "MIPGap" => 1E-9,
                                    "Presolve" => 1)
end

function get_juniper()
    return optimizer_with_attributes(Juniper.Optimizer, 
                                    MOI.Silent() => true, 
                                    "mip_solver" => get_cplex(), 
                                    "nl_solver" => get_ipopt())
end

function get_cplex()
    return optimizer_with_attributes(CPLEX.Optimizer, 
    MOI.Silent() => true, 
    "CPX_PARAM_PREIND" => 1,
    "CPX_PARAM_EPGAP" => 1E-9) 
end

function get_ipopt()
    return optimizer_with_attributes(Ipopt.Optimizer, 
                                    MOI.Silent() => true, 
                                    "sb" => "yes", 
                                    "max_iter" => Int(1E4))
end