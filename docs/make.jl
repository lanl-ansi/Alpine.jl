using Documenter, Apline

makedocs(
    format = :html,
    sitename = "Apline",
    pages = [
        "Introduction" => "intro.md",
        "How to Use" => "installation.md",
        "Choosing Sub-Solvers" => "choosingsolver.md",
        "Algorithm" => "algorithm.md",
        "Expression Guideline" => "expression.md",
        "Parameters" => "parameters.md",
        "Methods" => "functions.md",
        "Hacking Apline" => "hacking.md"
    ]
)


deploydocs(
    repo = "github.com/lanl-ansi/Apline.git",
)

