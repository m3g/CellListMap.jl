import Pkg; Pkg.add("Documenter")
using Documenter
import CellListMap
#DocMeta.setdocmeta!(CellListMap, :DocTestSetup, :(using CellListMap); recursive=true)
makedocs(
    modules=[CellListMap],
    sitename="CellListMap.jl",
    pages = [
        "Overview" => "index.md",
        "PeriodicSystems interface" => "PeriodicSystems.md",
        "Neighbor lists" => "neighborlists.md",
        "Low level interface" => "LowLevel.md",
        "Ecosystem integration" => "ecosystem.md",
        "From Python" => "python.md",
        "Help entries" => "help.md",
        "Citation" => "citation.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/CellListMap.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main", 
    versions = ["stable" => "v^", "v#.#" ],
)
