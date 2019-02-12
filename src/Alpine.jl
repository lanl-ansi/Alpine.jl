__precompile__()

module Alpine

using JuMP
using MathOptInterface
using Memento

# Create our module level logger (this will get precompiled)
const LOGGER = getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(Alpine)`
# NOTE: If this line is not included then the precompiled `Alpine.LOGGER` won't be registered at runtime.
__init__() = Memento.register(LOGGER)

"Suppresses information and warning messages output by PowerModels, for fine grained control use the Memento package"
function silence()
    info(LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    setlevel!(getlogger(Alpine), "error")
end


include("types.jl")
include("solver_options.jl")
include("MOI_wrapper.jl")

end # module
