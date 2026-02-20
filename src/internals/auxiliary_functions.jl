#=
    get_dim(unitcell, x, [y])

Try to determine the dimension (2 or 3) of the constructed system from the data
provided, error if imposible. 

=#
function get_dim(unitcell, x, y=nothing)
    Ds = ( 
        !isnothing(unitcell) ?  size(unitcell, 1) : 0,
        _element_length(x),
        _element_length(y),
    )
    D = maximum(Ds)
    if D == 0
        throw(ArgumentError("""\n
            Could not infer dimension (2 or 3) from the unitcell or coordinate arrays. 
            Got: unitcell: $(Ds[1]); x: $(Ds[2]); y: $(Ds[3])

        """))
    end
    if !(D in (2,3))
        throw(ArgumentError("""\n
            Dimension must be 2 or 3.
            Got: unitcell: $(Ds[1]); x: $(Ds[2]); y: $(Ds[3])

        """))
    end
    if any(d -> d != 0 && d != D, Ds)
        throw(ArgumentError("""\n
            Incompatible dimensions between unitcell and/or coordinates. 
            Got: unitcell: $(Ds[1]); x: $(Ds[2]); y: $(Ds[3])

        """))
    end
    return D
end
_element_length(::Nothing) = 0
_element_length(::AbstractVector{<:SVector{N}}) where {N} = N
_element_length(x::AbstractVector{<:AbstractVector}) = isempty(x) ? 0 : length(first(x))
_element_length(x::AbstractMatrix) = size(x,1)
_element_length(::AbstractVector) = throw(ArgumentError(" Coordinates must be vectors of vectors, or matrices."))

get_dim(sys::AbstractParticleSystem) = get_dim(sys.xpositions)
get_dim(::ParticleSystemPositions{N}) where {N} = N

#
# Auxiliary functions to control the exhibition of the progress meter
#
_next!(::Nothing) = nothing
_next!(p) = next!(p)

#
# Functions necessary for the projection/partition scheme
#
#=
    partition!(by, x::AbstractVector)

# Extended help

Function that reorders `x` vector by putting in the first positions the
elements with values satisfying `by(el)`. Returns the number of elements
that satisfy the condition.

=#
@inline function partition!(by, x::AbstractVector)
    iswap = 1
    @inbounds for i in eachindex(x)
        if by(x[i])
            if iswap != i
                x[iswap], x[i] = x[i], x[iswap]
            end
            iswap += 1
        end
    end
    return iswap - 1
end

#=
    project_particles!(projected_particles,cellⱼ,cellᵢ,Δc,Δc_norm,box)

# Extended help

Projects all particles of the cell `cellⱼ` into unnitary vector `Δc` with direction 
connecting the centers of `cellⱼ` and `cellᵢ`. Modifies `projected_particles`, and 
returns a view of `projected particles, where only the particles for which
the projection on the direction of the cell centers still allows the particle
to be within the cutoff distance of any point of the other cell.

=#
function project_particles!(
        projected_particles, cellⱼ, cellᵢ,
        Δc, Δc_norm, box::Box{UnitCellType, N}
    ) where {UnitCellType, N}
    if box.lcell == 1
        margin = box.cutoff + Δc_norm / 2 # half of the distance between centers
    else
        margin = box.cutoff * (1 + sqrt(N) / 2) # half of the diagonal of the cutoff-box
    end
    iproj = 0
    @inbounds for j in 1:cellⱼ.n_particles
        pⱼ = cellⱼ.particles[j]
        xproj = dot(pⱼ.coordinates - cellᵢ.center, Δc)
        if abs(xproj) <= margin
            iproj += 1
            projected_particles[iproj] = ProjectedParticle(pⱼ.index, xproj, pⱼ.coordinates, pⱼ.real)
        end
    end
    pp = @view(projected_particles[1:iproj])
    return pp
end
