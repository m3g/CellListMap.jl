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

export map_pairwise!, map_pairwise

# Testing file
const argon_pdb_file = joinpath("$(@__DIR__ )/../test/gromacs/argon/cubic.pdb")

# name holder  
function map_pairwise! end

"""
    map_pairwise(args...;kargs...) = map_pairwise!(args...;kargs...)

is an alias for `map_pairwise!` which is defined for two reasons: first, if the output of the funciton is immutable, it may be 
clearer to call this version, from a coding perspective. Second, the python interface through `juliacall` does not accept the 
bang as a valid character. 

"""
const map_pairwise = map_pairwise!

include("./linearalgebra.jl")
include("./show.jl")
include("./Box.jl")
include("./CellLists.jl")
include("./CellOperations.jl")
include("./ParticleSystem.jl")

# Core-computing infraestructure
include("./core_computing/auxiliary_functions.jl")
include("./core_computing/vicinal_cells.jl")
include("./core_computing/self.jl")
include("./core_computing/cross.jl")

# Utils
include("./neighborlists.jl")

#
# Test and example functions
#
include("./examples/examples.jl")
include("./testing.jl")

# Precompilation tools
include("precompile.jl")

end # module



