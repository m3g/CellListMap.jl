import Pkg; Pkg.add("Documenter")
using Documenter
using CellListMap
ENV["LINES"] = 10
ENV["COLUMNS"] = 120
makedocs(
    modules = [CellListMap],
    sitename = "CellListMap.jl",
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
        "Public Interface" => "API.md",
        #"From Python" => "python.md",
        "Citation" => "citation.md",
        "Migrating from 0.9" => "migrating.md",
        #"Internals" => "Internals.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/CellListMap.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#"],
)
