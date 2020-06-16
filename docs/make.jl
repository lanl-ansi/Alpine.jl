using Documenter, Alpine

makedocs(
    format = Documenter.HTML(analytics = "UA-367975-10", mathengine = Documenter.MathJax()),
    sitename = "Alpine",
    authors = "Harsha Nagarajan, Site Wang, Kaarthik Sundar, and contributors.",
    pages = [
        "Introduction" => "index.md",
        "How to Use" => "installation.md",
        "Choosing Sub-Solvers" => "choosingsolver.md",
        "Algorithm" => "algorithm.md",
        "Expression Guideline" => "expression.md",
        "Parameters" => "parameters.md",
        "Methods" => "functions.md",
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/Alpine.jl.git",
)

