"""
    NeighborPair{N,T,T2}

Structure that holds the information of a pair of particles that are neighbors within the cutoff distance.

## Fields accessed by the user:
- `i::Int`: index of the first particle in the original array of coordinates.
- `j::Int`: index of the second particle in the original array of coordinates.
- `x::SVector{N,T}`: coordinates of the first particle (minimum-image adjusted).
- `y::SVector{N,T}`: coordinates of the second particle (minimum-image adjusted).
- `d::T`: Euclidean distance between the particles (computed lazily).
- `d2::T2`: squared Euclidean distance between the particles.

"""
struct NeighborPair{N, T, T2}
    i::Int
    j::Int
    x::SVector{N, T}
    y::SVector{N, T}
    d2::T2
end

import Base: getproperty, propertynames
getproperty(p::NeighborPair, s::Symbol) = getproperty(p, Val(s))
getproperty(p::NeighborPair, ::Val{S}) where {S} = getfield(p, S)
# Lazy computation of the distance
getproperty(p::NeighborPair, ::Val{:d}) = sqrt(getfield(p, :d2))
# Property names
propertynames(::NeighborPair) = (:i, :j, :x, :y, :d, :d2)
