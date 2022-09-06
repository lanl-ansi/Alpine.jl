#---------------------------;
# MIP solvers - commercial  ;
#---------------------------;
# https://github.com/jump-dev/Gurobi.jl
function get_gurobi()
    GRB_ENV = Gurobi.Env()
    return optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), # To avoid printing License info multiple times
        MOI.Silent() => true,
        "Presolve" => 1,
    )
end

# https://github.com/jump-dev/CPLEX.jl
function get_cplex()
    return optimizer_with_attributes(
        CPLEX.Optimizer,
        MOI.Silent() => true,
        "CPX_PARAM_PREIND" => 1,
    )
end

#---------------------------;
# MIP solvers - open-source
#---------------------------;
# https://github.com/jump-dev/HiGHS.jl
function get_highs()
    return JuMP.optimizer_with_attributes(
        HiGHS.Optimizer,
        "presolve" => "on",
        "log_to_console" => false,
    )
end

# https://github.com/jump-dev/Cbc.jl
function get_cbc()
    return optimizer_with_attributes(Cbc.Optimizer, MOI.Silent() => true)
end

#---------------------------;
# Continuous local solver   ; 
#---------------------------;
# https://github.com/jump-dev/Ipopt.jl
function get_ipopt()
    return optimizer_with_attributes(
        Ipopt.Optimizer,
        MOI.Silent() => true,
        "sb" => "yes",
        "max_iter" => Int(1E4),
    )
end

#----------------------;
# Convex MINLP solver  ;
#----------------------;
# https://github.com/jump-dev/Pavito.jl
function get_pavito(mip_solver, cont_solver)
    return optimizer_with_attributes(
        Pavito.Optimizer,
        MOI.Silent() => true,
        "mip_solver" => mip_solver,
        "cont_solver" => cont_solver,
        "mip_solver_drives" => false,
    )
end

#-------------------------------;
# Non-convex Local MINLP solver ;
#-------------------------------;
# https://github.com/lanl-ansi/Juniper.jl
function get_juniper(mip_solver, nl_solver)
    return optimizer_with_attributes(
        Juniper.Optimizer,
        MOI.Silent() => true,
        "mip_solver" => mip_solver,
        "nl_solver" => nl_solver,
    )
end
