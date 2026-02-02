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

export NeighborPair
export map_pairwise!, map_pairwise

# Testing file
const argon_pdb_file = joinpath("$(@__DIR__ )/../test/gromacs/argon/cubic.pdb")

"""
    NeighborPair{N,T}

Structure that holds the information of a pair of particles that are neighbors
within the cutoff distance.

## Fields accessed by the user:
- `i::Int`: index of the first particle in the original array of coordinates.
- `j::Int`: index of the second particle in the original array of coordinates.
- `x::SVector{N,T}`: coordinates of the first particle (minimum-image adjusted).
- `y::SVector{N,T}`: coordinates of the second particle (minimum-image adjusted).
- `d::T`: Euclidean distance between the particles (computed lazily).
- `d2::T`: squared Euclidean distance between the particles.

## Example

```julia-repl
julia> sys = ParticleSystem(positions=rand(SVector{3,Float64},100), cutoff=0.1, unitcell=[1,1,1], output=0.0);

julia> map_pairwise!((pair, out) -> out += 1/pair.d, sys)
```

"""
struct NeighborPair{N,T,T2}
    i::Int
    j::Int
    x::SVector{N,T}
    y::SVector{N,T}
    d2::T2
end
Base.getproperty(p::NeighborPair, s::Symbol) = getproperty(p, Val(s))
Base.getproperty(p::NeighborPair, ::Val{S}) where {S} = getfield(p, S)
# Lazy computation of the distance
Base.getproperty(p::NeighborPair, ::Val{:d}) = sqrt(getfield(p, :d2))
# Property names
Base.propertynames(::NeighborPair) = (:i, :j, :x, :y, :d, :d2)

# name holder
function map_pairwise! end

"""
    map_pairwise!(args...;kargs...) = map_pairwise!(args...;kargs...)

is an alias for `map_pairwise!` which is defined for two reasons: first, if the output of the function is immutable, it may be
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



