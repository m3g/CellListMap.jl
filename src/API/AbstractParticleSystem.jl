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

Thin wrapper of particle coordinates, containing a `Vector{SVector{N,T}}` that carries the coordinates
of the system and a `updated` flag, to instruct the update, or not, of the cell lists. The fields
of `ParticleSystemPositions` are internal. Coordinates can be fetched by `getindex` and updated
with `setindex!`, as usual arrays.

ParticleSystemPositions can be initialized with vectors of vectors of coordinates, or with (M,N) matrices, 
where M is the number of particles and N is the dimension of the system.

"""
struct ParticleSystemPositions{N,T,V<:AbstractVector{SVector{N,T}}}
    x::V
    updated::Ref{Bool}
end

function ParticleSystemPositions(x::AbstractVector{<:AbstractVector})
    M = length(x)
    N, T = if M > 0
        length(first(x)), eltype(first(x))
    else
        # This will error for empty vectors of non-fixed-size element vectors
        length(eltype(x)), eltype(eltype(x))
    end
    x_static = Vector{SVector{N,T}}(undef, M)
    for i in eachindex(x, x_static)
        x_static[i] = SVector{N,T}(ntuple(j -> x[i][j], N))
    end
    ParticleSystemPositions{N,T,typeof(x_static)}(x_static, Ref(true))
end

# If the the elements are ot fixed-sized and the vector is empty, pass DIM
function ParticleSystemPositions(x::AbstractVector{<:AbstractVector}, DIM)
    M = length(x)
    T = eltype(eltype(x))
    x_static = Vector{SVector{DIM,T}}(undef, M)
    for i in eachindex(x, x_static)
        x_static[i] = SVector{DIM,T}(ntuple(j -> x[i][j], DIM))
    end
    ParticleSystemPositions{DIM,T,typeof(x_static)}(x_static, Ref(true))
end

function ParticleSystemPositions(x::AbstractMatrix{T}, DIM=size(x,1)) where {T}
    M = size(x,2)
    xv = Vector{SVector{DIM,T}}(undef, M)
    for i in eachindex(xv)
        xv[i] = SVector{DIM,T}(@view(x[:,i]))
    end
    ParticleSystemPositions{DIM,T,typeof(xv)}(xv, Ref(true))
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
Base.append!(p::ParticleSystemPositions, x) = (p.updated[] = true; append!(p.x, x); p)
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
