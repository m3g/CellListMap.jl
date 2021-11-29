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
        "Units, autodiff, etc." => "units_etc.md",
        "Periodic conditions" => "pbc.md",
        "Parallelization" => "parallelization.md",
        "Performance" => "performance.md",
        "Options" => "options.md",
        "From Python" => "python.md",
        "Help entries" => "help.md",
        "Reference" => "reference.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/CellListMap.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main", 
    versions = ["stable" => "v^", "v#.#" ],
)
