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

$(TYPEDEF)

Structure that carries the information necessary for `pairwise!` computations,
for systems with one set of positions (thus, replacing the loops over `N(N-1)` 
pairs of particles of the set). Can be used to control dispatch.

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

$(TYPEDEF)

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
