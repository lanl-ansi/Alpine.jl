using Documenter, Alpine
using DocumenterTools: Themes

for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "src/assets/themes/style.scss"), String)
    theme = read(joinpath(@__DIR__, "src/assets/themes/$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "src/assets/themes/$(w).scss"), header * "\n" * theme)
end

Themes.compile(
    joinpath(@__DIR__, "src/assets/themes/documenter-light.css"),
    joinpath(@__DIR__, "src/assets/themes/light.scss"),
)
Themes.compile(
    joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"),
    joinpath(@__DIR__, "src/assets/themes/dark.scss"),
)

makedocs(
    sitename = "Alpine",
    authors = "Harsha Nagarajan, Site Wang, Kaarthik Sundar and contributors",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    # strict = true,
    pages = [
        "Introduction" => "index.md",
        "Choosing Sub-solvers" => "choosingsolver.md",
        # "Algorithm" => "algorithm.md",
        "Solver Options" => "parameters.md",
        "Methods" => "functions.md",
    ],
)

deploydocs(repo = "github.com/lanl-ansi/Alpine.jl.git", push_preview = true)
