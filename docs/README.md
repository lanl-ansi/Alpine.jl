# Documentation for Alpine

## Installation
Alpine.jl's documentation is generated using [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). To install it, run the following command in a Julia session:

```julia
Pkg.add("Documenter")
Pkg.add("DocumenterTools")
```

## Building the Docs
To preview the `html` output of the documents, run the following command:

```julia
julia --color=yes make.jl
```

You can then view the documents in `build/index.html`.