using Documenter, POD

makedocs(
    format = :html,
    sitename = "POD",
    pages = [
        "Home" => "index.md",
        "Getting Started" =>
            [
            "Installation" => "installation.md",
            "Basic Usage" => "usage.md",
            "Choosing Sub-Solvers" => "choosingsolver.md"
            ],
        "Parameters" => "parameters.md",
        "Expression Guideline" => "expression.md",
        "Hacking POD" => "hacking.md"
    ]
)


deploydocs(
    repo = "github.com/lanl-ansi/POD.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
