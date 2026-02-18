"""
    AbstractParticleSystem{OutputName}

An abstract type representing a particle system that computes interactions
between particles using cell lists. Can be used to control dispatch of functions
that operate on particle systems, because all particle systems should be a subtype
of this type.

Use the `ParticleSystem` constructor and interface for anything else than dispatch.

"""
abstract type AbstractParticleSystem{OutputName} end

"""
    ParticleSystemPositions{N,T}

Wrapper around a `Vector{SVector{N,T}}` that carries particle coordinates and an internal
`updated` flag. When coordinates are mutated through the supported interface, the flag is
set automatically, so that cell lists are recomputed on the next call to `pairwise!`.

A `ParticleSystemPositions` can be constructed from:
- A vector of vectors (e.g. `Vector{Vector{Float64}}`).
- A vector of `SVector`s (e.g. `Vector{SVector{3,Float64}}`).
- An `(D, M)` matrix, where `D` is the dimension and `M` the number of particles.

## Mutating interface

The following functions mutate the positions **and** flag the array as updated, triggering
recomputation of the cell lists on the next `pairwise!` call:

| Function      | Description                                     |
|:------------- |:----------------------------------------------- |
| `setindex!`   | Set the position of a single particle by index   |
| `empty!`      | Remove all positions                             |
| `resize!`     | Resize the number of positions                   |
| `append!`     | Append positions from another collection         |
| `push!`       | Add element                                      |
| Broadcasting  | In-place broadcast (e.g. `p .= new_positions`)   |

## Read-only interface

| Function      | Description                                      |
|:------------- |:------------------------------------------------ |
| `getindex`    | Retrieve the position of a particle by index      |
| `length`      | Number of particles                               |
| `size`        | Size tuple `(length,)`                            |
| `axes`        | Index axes of the underlying vector               |
| `keys`        | Linear indices                                    |
| `eachindex`   | Iterator over valid indices                       |
| `firstindex`  | First valid index                                 |
| `lastindex`   | Last valid index                                  |
| `first`       | First position                                    |
| `last`        | Last position                                     |
| `ndims`       | Always returns `1`                                |
| `iterate`     | Iteration protocol                                |
| `copy`        | Shallow copy (preserves `updated` flag)           |
| `similar`     | Allocate an uninitialized array of same shape      |
| `view`        | Create a view sharing the `updated` flag          |
| `show`        | Pretty-printing                                   |

"""
struct ParticleSystemPositions{N,T,V<:AbstractVector{SVector{N,T}}}
    x::V
    updated::Ref{Bool}
end

# Single-arg constructors: infer DIM and dispatch through Val for type stability
function ParticleSystemPositions(x::AbstractVector{SVector{N,T}}) where {N,T}
    x_copy = copy(x)
    ParticleSystemPositions{N,T,typeof(x_copy)}(x_copy, Ref(true))
end

function ParticleSystemPositions(x::AbstractVector{<:AbstractVector})
    N = length(first(x))
    ParticleSystemPositions(x, Val(N))
end

function ParticleSystemPositions(x::AbstractMatrix)
    ParticleSystemPositions(x, Val(size(x, 1)))
end

# Two-arg with Val{N}: type-stable constructors
function ParticleSystemPositions(x::AbstractVector{SVector{N,T}}, ::Val) where {N,T}
    x_copy = copy(x)
    ParticleSystemPositions{N,T,typeof(x_copy)}(x_copy, Ref(true))
end

function ParticleSystemPositions(x::AbstractVector{<:AbstractVector}, ::Val{N}) where {N}
    T = eltype(eltype(x))
    M = length(x)
    x_static = Vector{SVector{N,T}}(undef, M)
    for i in eachindex(x, x_static)
        x_static[i] = SVector{N,T}(ntuple(j -> x[i][j], Val(N)))
    end
    ParticleSystemPositions{N,T,typeof(x_static)}(x_static, Ref(true))
end

function ParticleSystemPositions(x::AbstractMatrix{T}, ::Val{N}) where {T,N}
    M = size(x, 2)
    xv = Vector{SVector{N,T}}(undef, M)
    for i in eachindex(xv)
        xv[i] = SVector{N,T}(@view(x[:, i]))
    end
    ParticleSystemPositions{N,T,typeof(xv)}(xv, Ref(true))
end

function Base.empty!(p::ParticleSystemPositions)
    p.updated[] = true
    empty!(p.x)
    return p
end
function Base.resize!(p::ParticleSystemPositions, n::Integer)
    p.updated[] = true
    resize!(p.x, n)
    return p
end
function Base.push!(p::ParticleSystemPositions{N,T}, v::SVector{N,T}) where {N,T}
    p.updated[] = true
    push!(p.x, v)
    return p
end

Base.getindex(p::ParticleSystemPositions, i) = p.x[i]
function Base.setindex!(p::ParticleSystemPositions{N,T}, v::SVector{N,T}, i) where {N,T}
    p.updated[] = true
    p.x[i] = v
end
Base.length(p::ParticleSystemPositions) = length(p.x)
Base.ndims(p::ParticleSystemPositions) = 1
Base.ndims(::Type{<:ParticleSystemPositions}) = 1
Base.iterate(p::ParticleSystemPositions, i=firstindex(p.x)) = i > length(p.x) ? nothing : (p.x[i], i + 1)
Base.axes(p::ParticleSystemPositions) = axes(p.x)
Base.keys(p::ParticleSystemPositions) = LinearIndices(p.x)
Base.size(p::ParticleSystemPositions) = size(p.x)
Base.firstindex(p::ParticleSystemPositions) = firstindex(p.x)
Base.lastindex(p::ParticleSystemPositions) = lastindex(p.x)
Base.first(p::ParticleSystemPositions) = first(p.x)
Base.last(p::ParticleSystemPositions) = last(p.x)
Base.copy(p::ParticleSystemPositions{N,T,V}) where {N,T,V<:Vector} = ParticleSystemPositions{N,T,V}(copy(p.x), Ref(p.updated[]))
function Base.append!(p::ParticleSystemPositions, x) 
    p.updated[] = true 
    append!(p.x, x)
    return p
end
Base.similar(p::ParticleSystemPositions{N,T}) where {N,T} = ParticleSystemPositions{N,T,Vector{SVector{N,T}}}(similar(p.x), Ref(true))
Base.eachindex(p::ParticleSystemPositions) = eachindex(p.x)

# Broadcast interface
struct BroadcastParticleSystemPositions <: Broadcast.BroadcastStyle end
Base.BroadcastStyle(::Type{<:ParticleSystemPositions}) = BroadcastParticleSystemPositions()
Base.BroadcastStyle(::BroadcastParticleSystemPositions, ::Broadcast.DefaultArrayStyle{0}) = BroadcastParticleSystemPositions()
Base.BroadcastStyle(::BroadcastParticleSystemPositions, ::Broadcast.DefaultArrayStyle{1}) = BroadcastParticleSystemPositions()
Base.Broadcast.broadcastable(p::ParticleSystemPositions) = p

# Out-of-place broadcast returns a plain Array
function Base.similar(bc::Broadcast.Broadcasted{BroadcastParticleSystemPositions}, ::Type{T}, axes) where {T}
    similar(Array{T}, axes)
end

# In-place broadcast into ParticleSystemPositions
@inline function Base.copyto!(dest::ParticleSystemPositions, bc::Broadcast.Broadcasted)
    dest.updated[] = true
    @inbounds @simd for i in eachindex(dest.x)
        dest.x[i] = bc[i]
    end
    return dest
end
function Base.copyto!(dest::ParticleSystemPositions, src::ParticleSystemPositions)
    dest.updated[] = true
    copyto!(dest.x, src.x)
    return dest
end

function Base.view(p::ParticleSystemPositions{N,T}, inds...) where {N,T}
    x_view = view(p.x, inds...)
    ParticleSystemPositions{N,T,typeof(x_view)}(x_view, p.updated[])
end
function Base.show(io::IO, ::MIME"text/plain", p::ParticleSystemPositions) 
    print(io, "ParticleSystemPositions, ")
    print(io, "updated: ", p.updated[], ", with ")
    show(IOContext(io, :title => false), MIME"text/plain"(), p.x)
end

"""

    ParticleSystem1

Structure that carries the information necessary for `pairwise!` computations,
for systems with one set of positions (thus, replacing the loops over `N(N-1)` 
pairs of particles of the set). Can be used to control dispatch.

Use the `ParticleSystem` constructor and interface for anything else than dispatch.

"""
mutable struct ParticleSystem1{OutputName, V, O, B, C, A, VC} <: AbstractParticleSystem{OutputName}
    xpositions::V
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
    validate_coordinates::VC
end

"""

    ParticleSystem2

Structure that carries the information necessary for `pairwise!` computations,
for systems with two set of positions (thus, replacing the loops over `N×M` 
pairs of particles, being `N` and `M` the number of particles of each set).
Can be used to control dispatch.

Use the `ParticleSystem` constructor and interface for anything else than dispatch.

"""
mutable struct ParticleSystem2{OutputName, V, O, B, C, A, VC} <: AbstractParticleSystem{OutputName}
    xpositions::V
    ypositions::V
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
    validate_coordinates::VC
end
ParticleSystem2{OutputName}(vx::V, vy::V, o::O, b::B, c::C, vo::Vector{O}, a::A, p::Bool, vc::VC) where {OutputName, V, O, B, C, A, VC} =
    ParticleSystem2{OutputName, V, O, B, C, A, VC}(vx, vy, o, b, c, vo, a, p, vc)
