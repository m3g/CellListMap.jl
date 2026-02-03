module CellListMap

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using TestItems: @testitem, @testmodule
using Compat: @compat
using ProgressMeter: Progress, next!
using Parameters: @unpack, @with_kw
using StaticArrays: SVector, SMatrix, @SVector, @SMatrix, MVector, MMatrix, FieldVector
using Setfield: @set!
using LinearAlgebra: cross, diagm, I
using Base.Threads: nthreads, @spawn
using Base: @lock # not exported in 1.6
using ChunkSplitters: index_chunks, RoundRobin, Consecutive

# Exported names of ParticleSystem interface
export ParticleSystem
export NeighborPair
export pairwise!
export update_cutoff!
export update_unitcell!
export resize_output!
@compat public copy_output, reset_output!, reset_output, reducer, reducer!
@compat public wrap_relative_to

# Specific for neighborlist interface
export InPlaceNeighborList
export update!
export neighborlist, neighborlist!

include("./API/neighborpair.jl")
include("./internals/linearalgebra.jl")
include("./internals/show.jl")
include("./internals/Box.jl")
include("./internals/CellLists.jl")
include("./internals/CellOperations.jl")
include("./API/ParticleSystem.jl")

# Core-computing infraestructure
include("./internals/auxiliary_functions.jl")
include("./internals/vicinal_cells.jl")
include("./internals/self.jl")
include("./internals/cross.jl")

# Neighborlists
include("./API/neighborlists.jl")

# Precompilation tools
include("./internals/precompile.jl")

# Test file used in some doc strings and example blocs
const argon_pdb_file = joinpath("$(@__DIR__ )/../test/applications/gromacs/argon/cubic.pdb")

end # module



