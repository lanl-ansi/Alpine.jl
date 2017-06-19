using Documenter, POD

makedocs(
    format = :html,
    sitename = "POD",
    pages = [
        "Introduction" => "intro.md",
        "Installation" => "installation.md",
        "Demo" => "demo.md",
        "Parameters" => "parameters.md",
        "Choosing Solvers" => "choosingsolver.md",
        "Methods" => "functions.md"
    ]
)


deploydocs(
    repo = "github.com/lanl-ansi/POD.git",
    target = "build",
    branch = "monomial",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
