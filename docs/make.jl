using Documenter, ControlSystems, Plots

makedocs(modules=[POD])

function myDeps()
    return
end

deploydocs(
    repo = "github.com/lanl-ansi/POD.git",
    latest = "master",
    julia = "0.5",
    deps = myDeps
)
