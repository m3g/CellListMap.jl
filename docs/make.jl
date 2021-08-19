import Pkg
Pkg.add("Documenter")
using Documenter
using CellListMap
DocMeta.setdocmeta!(CellListMap, :DocTestSetup, :(using CellListMap); recursive=true)
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[CellListMap],
    sitename="CellListMap.jl",
    pages = [
        "Overview" => "index.md",
        "Examples" => "examples.md",
        "Periodic conditions" => "pbc.md",
        "Parallelization" => "parallelization.md",
        "Performance" => "performance.md",
        "Reference" => "reference.md",
        "Help entries" => "help.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/CellListMap.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
