module CellListMap

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using TestItems: @testitem, @testmodule
using Compat: @compat
using ProgressMeter: Progress, next!
using StaticArrays: SVector, SMatrix, @SVector, @SMatrix, MVector, MMatrix, FieldVector
using Setfield: @set!
using LinearAlgebra: cross, diagm, I, dot, norm
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
@compat public AbstractParticleSystem, ParticleSystem1, ParticleSystem2
@compat public copy_output, reset_output!, reset_output, reducer, reducer!, reduce_output!
@compat public wrap_relative_to, get_computing_box

# Specific for neighborlist interface
export InPlaceNeighborList
export update!
export neighborlist, neighborlist!

# NeighborPair structure
include("./API/NeighborPair.jl")

# ParticleSystemPositions type (needed by CellLists.jl)
include("./API/AbstractParticleSystem.jl")
include("./API/ParticleSystemPositions.jl")

# Core-computing
include("./internals/show.jl")
include("./internals/Box.jl")
include("./internals/CellLists.jl")
include("./internals/CellOperations.jl")
include("./internals/auxiliary_functions.jl")
include("./internals/vicinal_cells.jl")
include("./internals/self.jl")
include("./internals/cross.jl")
include("./internals/ParticleSystem.jl")
include("./API/ParticleSystem.jl")
include("./API/parallel_custom.jl")
include("./API/updating.jl")
include("./API/pairwise.jl")
include("./API/get_computing_box.jl")

# Neighborlists
include("./internals/neighborlist.jl")
include("./API/neighborlist.jl")

# Precompilation tools
include("./internals/precompile.jl")

# Test file used in some doc strings and example blocs
const argon_pdb_file = joinpath("$(@__DIR__)/../test/applications/gromacs/argon/cubic.pdb")

end # module
