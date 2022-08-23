# MIP solver - commercial 
function get_gurobi()
    return optimizer_with_attributes(
        Gurobi.Optimizer,
        MOI.Silent() => true,
        "Presolve" => 1,
    )
end

function get_cplex()
    return optimizer_with_attributes(
        CPLEX.Optimizer,
        MOI.Silent() => true,
        "CPX_PARAM_PREIND" => 1,
    )
end

# MIP solver - open-source
function get_cbc()
    return optimizer_with_attributes(Cbc.Optimizer, MOI.Silent() => true)
end

# Local solver
function get_ipopt()
    return optimizer_with_attributes(
        Ipopt.Optimizer,
        MOI.Silent() => true,
        "sb" => "yes",
        "max_iter" => Int(1E4),
    )
end

# Convex MINLP solver
function get_pavito(mip_solver, cont_solver)
    return optimizer_with_attributes(
        Pavito.Optimizer,
        MOI.Silent() => true,
        "mip_solver" => mip_solver,
        "cont_solver" => cont_solver,
        "mip_solver_drives" => false,
    )
end

# Non-convex Local MINLP solver
function get_juniper(mip_solver, nl_solver)
    return optimizer_with_attributes(
        Juniper.Optimizer,
        MOI.Silent() => true,
        "mip_solver" => mip_solver,
        "nl_solver" => nl_solver,
    )
end
