using Documenter, POD

makedocs(
    format = :html,
    sitename = "POD",
    pages = [
        "Introduction" => "intro.md",
        "How to Use" => "installation.md",
        "Choosing Sub-Solvers" => "choosingsolver.md",
        "Algorithm" => "algorithm.md",
        "Expression Guideline" => "expression.md",
        "Parameters" => "parameters.md",
        "Methods" => "functions.md",
        "Hacking POD" => "hacking.md"
    ]
)


deploydocs(
    repo = "github.com/lanl-ansi/POD.git",
)

