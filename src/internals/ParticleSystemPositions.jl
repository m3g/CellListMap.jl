#=
    ParticleSystemPositions{N,T}

Internal type. Wrapper around a `Vector{SVector{N,T}}` that carries particle coordinates
and an `updated::Ref{Bool}` flag. Any mutation through the standard array interface
(indexing, `push!`, `resize!`, `empty!`, `append!`, in-place broadcast) sets the flag,
triggering cell-list recomputation on the next `pairwise!` call.

Users interact with this type only through `system.xpositions` / `system.ypositions`,
treating it as an ordinary vector. To replace the entire coordinate array or update
other system properties, use `update!`.

=#
struct ParticleSystemPositions{N,T,V<:AbstractVector{SVector{N,T}}}
    x::V
    updated::Ref{Bool}
end

function ParticleSystemPositions(x::AbstractVector{SVector{N,T}}) where {N,T}
    x_copy = copy(x)
    ParticleSystemPositions{N,T,typeof(x_copy)}(x_copy, Ref(true))
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
    ParticleSystemPositions{N,T,typeof(x_view)}(x_view, p.updated)
end

function Base.show(io::IO, ::MIME"text/plain", p::ParticleSystemPositions) 
    print(io, "ParticleSystemPositions, ")
    print(io, "updated: ", p.updated[], ", with ")
    show(IOContext(io, :title => false), MIME"text/plain"(), p.x)
end
