import Pkg; Pkg.add("Documenter")
using Documenter
using CellListMap
#DocMeta.setdocmeta!(CellListMap, :DocTestSetup, :(
#    #using CellListMap
#    Base.active_repl.options.iocontext[:displaysize] = (9, 80)
#); recursive=true)
makedocs(
    modules=[CellListMap],
    sitename="CellListMap.jl",
    pages = [
        "Overview" => "index.md",
        "Neighbor lists" => "neighborlists.md",
        "ParticleSystem" => "ParticleSystem/introduction.md",
        "Single set: Simple outputs" => "ParticleSystem/single_set_simple.md",
        "Single set: Compound outputs" => "ParticleSystem/single_set_compound.md",
        "Two sets of particles" => "ParticleSystem/two_sets.md",
        "Updating the system" => "ParticleSystem/updating.md",
        "Options" => "ParticleSystem/options.md",
        "Complete examples" => "ParticleSystem/examples.md",
        "Unitcell requirements" => "unitcell.md",
        "Ecosystem integration" => "ecosystem.md",
        "From Python" => "python.md",
        "Citation" => "citation.md",
        "Internals" => "Internals.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/CellListMap.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main", 
    versions = ["stable" => "v^", "v#.#" ],
)
