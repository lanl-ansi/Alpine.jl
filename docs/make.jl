using Documenter, Alpine

makedocs(
    format = :html,
    sitename = "Alpine",
    pages = [
        "Introduction" => "intro.md",
        "How to Use" => "installation.md",
        "Choosing Sub-Solvers" => "choosingsolver.md",
        "Algorithm" => "algorithm.md",
        "Expression Guideline" => "expression.md",
        "Parameters" => "parameters.md",
        "Methods" => "functions.md",
        "Hacking Alpine" => "hacking.md"
    ]
)


deploydocs(
    repo = "github.com/lanl-ansi/Alpine.jl.git",
)

