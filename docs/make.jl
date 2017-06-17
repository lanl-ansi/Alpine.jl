using Documenter, POD

makedocs(
    format = :html,
    sitename = "POD",
    pages = [
        "Introduction" => "intro.md",
        "Installation" => "installation.md",
        "Demo" => "firstdemo.md",
        "Choosing Solvers" => "choosingsolver.md"
    ]
)


deploydocs(
    repo = "github.com/lanl-ansi/POD.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
