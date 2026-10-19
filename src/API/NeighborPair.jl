"""
    NeighborPair{N,T,T2}

Structure that holds the information of a pair of particles that are neighbors within the cutoff distance.

## Fields accessed by the user:
- `i::Int`: index of the first particle in the original array of coordinates.
- `j::Int`: index of the second particle in the original array of coordinates.
- `x::SVector{N,T}`: coordinates of the first particle (wrapped into the unit cell).
- `y::SVector{N,T}`: coordinates of the image of the second particle within the cutoff of `x`.
- `d::T`: Euclidean distance between the particles (computed lazily).
- `d2::T2`: squared Euclidean distance between the particles.
- `shift::SVector{N,Int}`: lattice translation of particle `j` relative to particle `i`,
  in units of the unit cell vectors (see below).

Specifically, `pair.d == norm(pair.y - pair.x)`, and `pair.d` will be always equal or 
smaller than the cutoff, as only pairs within the cutoff are considered. `pair.y`
is the position of the image of particle `pair.j` within the cutoff of `pair.x`. If the 
unit cell is large relative to the cutoff, this is the minimum image; otherwise, the 
function is called once for each image within the cutoff, possibly including
images of particle `i` itself (with `pair.i == pair.j`).

The `shift` field gives the periodic image of particle `j` that is within the cutoff
of particle `i`, relative to the original (input) coordinates of the particles. That is,
if `x` and `y` are the input coordinate arrays (`y = x` for single-set computations),
and `unitcell` is the unit cell matrix (or `Diagonal(sides)` for orthorhombic cells), then

    y[pair.j] + unitcell * pair.shift - x[pair.i] ≈ pair.y - pair.x

For non-periodic systems the shift is always zero.

"""
struct NeighborPair{N, T, T2}
    i::Int
    j::Int
    x::SVector{N, T}
    y::SVector{N, T}
    d2::T2
    shift::SVector{N, Int}
end
NeighborPair(i, j, x::SVector{N, T}, y::SVector{N, T}, d2) where {N, T} =
    NeighborPair(i, j, x, y, d2, zero(SVector{N, Int}))

# Lattice translation of particle `pⱼ` relative to particle `pᵢ` (see `ParticleWithIndex`).
@inline _pair_shift(pᵢ, pⱼ) = SVector{length(pⱼ.shift), Int}(pⱼ.shift) - SVector{length(pᵢ.shift), Int}(pᵢ.shift)

# In self-computations where all neighboring cells are visited (triclinic cells), each pair
# is found from both particles, and only the one with `i < j` is kept. A particle and its own
# periodic image (possible when the cell is small relative to the cutoff) are found
# with opposite shifts, and only the one with a lexicographically positive shift is kept
# (which also skips the particle with itself, with a zero shift).
@inline function _skip_self_pair(pᵢ, pⱼ)
    pᵢ.index < pⱼ.index && return false
    pᵢ.index > pⱼ.index && return true
    return !_lexicographically_positive(pⱼ.shift - pᵢ.shift)
end
@inline function _lexicographically_positive(s::SVector{N}) where {N}
    for k in 1:N
        s[k] > 0 && return true
        s[k] < 0 && return false
    end
    return false
end

import Base: getproperty, propertynames
getproperty(p::NeighborPair, s::Symbol) = getproperty(p, Val(s))
getproperty(p::NeighborPair, ::Val{S}) where {S} = getfield(p, S)
# Lazy computation of the distance
getproperty(p::NeighborPair, ::Val{:d}) = sqrt(getfield(p, :d2))
# Property names
propertynames(::NeighborPair) = (:i, :j, :x, :y, :d, :d2, :shift)
