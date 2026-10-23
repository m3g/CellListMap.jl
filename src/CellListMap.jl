module CellListMap

using TestItems: @testitem, @testmodule
using Compat: @compat
using ProgressMeter: Progress, next!
using StaticArrays: SVector, SMatrix, @SVector, @SMatrix, MVector, MMatrix, FieldVector, setindex
using Setfield: @set!
using LinearAlgebra: cross, diagm, I, dot, norm
using Base.Threads: nthreads, @spawn, Atomic, atomic_add!
using Base: @lock # not exported in 1.6
using ChunkSplitters: index_chunks, Consecutive

# Exported names of ParticleSystem interface
export ParticleSystem
export NeighborPair
export pairwise!
export resize_output!
@compat public AbstractParticleSystem, ParticleSystem1, ParticleSystem2
@compat public ParticleSystemPositions
@compat public copy_output, reset_output!, reset_output, reducer, reducer!, reduce_output!
@compat public wrap_relative_to, get_computing_box, get_nbatches

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
include("./internals/pairs.jl")
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

# test_file(name): path of the test file `name` (e.g. "cubic.pdb", the argon system used in
# the examples). The test files are a lazy artifact (see Artifacts.toml), downloaded on the
# first call. The method is defined in the TestFilesExt extension, loaded with
# `using LazyArtifacts`.
function test_file end

end # module
