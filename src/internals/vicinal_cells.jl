#
# Compute interactions between vicinal cells
#
function _vicinal_cell_interactions!(f::F, box::Box, cellᵢ, cellⱼ, cl::CellList{N, T}, output, ibatch, ::Val{Skip} = Val(false)) where {F <: Function, N, T, Skip}
    # project particles in vector connecting cell centers
    Δc = cellⱼ.center - cellᵢ.center
    Δc_norm = norm(Δc)
    Δc = Δc / Δc_norm
    pp = project_particles!(cl.projected_particles[ibatch], cellⱼ, cellᵢ, Δc, Δc_norm, box)
    if length(pp) > 0
        output = _vinicial_cells!(f, box, cellᵢ, pp, Δc, output, Val(Skip))
    end
    return output
end

#
# The criteria form skipping computations is different then in Orthorhombic or Triclinic boxes. The
# first parameter (the type Self, or Cross, computation, is not needed here, because the symmetry
# allows to never compute repeated interactions anyway.
#
function _vinicial_cells!(f::F, box::Box{<:OrthorhombicCellType}, cellᵢ, pp, Δc, output, ::Val{Skip}) where {F <: Function, Skip}
    (; cutoff, cutoff_sqr, inv_rotation) = box
    # Loop over particles of cell icell
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cellᵢ.particles[i]
        # project particle in vector connecting cell centers
        xpᵢ = pᵢ.coordinates
        xpᵢ_rot = inv_rotation * xpᵢ
        xproj = dot(xpᵢ - cellᵢ.center, Δc)
        # Partition pp array according to the current projections
        n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)
        # Compute the interactions
        for j in 1:n
            @inbounds pⱼ = pp[j]
            (pᵢ.real | pⱼ.real) || continue
            xpⱼ = pⱼ.coordinates
            d2 = sum(abs2, xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                pair = NeighborPair(pᵢ.index, pⱼ.index, xpᵢ_rot, inv_rotation * xpⱼ, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end

# Here Skip determines if the interactions are self or cross, in such a way
# that, for self-computations, we need to skip the interactions when i >= j.
function _vinicial_cells!(f::F, box::Box{<:TriclinicCell}, cellᵢ, pp, Δc, output, ::Val{Skip}) where {F <: Function, Skip}
    (; cutoff, cutoff_sqr, inv_rotation) = box
    # Loop over particles of cell icell
    for i in 1:cellᵢ.n_particles
        @inbounds pᵢ = cellᵢ.particles[i]
        pᵢ.real || continue
        # project particle in vector connecting cell centers
        xpᵢ = pᵢ.coordinates
        xpᵢ_rot = inv_rotation * xpᵢ
        xproj = dot(xpᵢ - cellᵢ.center, Δc)
        # Partition pp array according to the current projections
        n = partition!(el -> abs(el.xproj - xproj) <= cutoff, pp)
        # Compute the interactions
        for j in 1:n
            @inbounds pⱼ = pp[j]
            if Skip
                pᵢ.index >= pⱼ.index && continue
            end
            xpⱼ = pⱼ.coordinates
            d2 = sum(abs2, xpᵢ - xpⱼ)
            if d2 <= cutoff_sqr
                pair = NeighborPair(pᵢ.index, pⱼ.index, xpᵢ_rot, inv_rotation * xpⱼ, d2)
                output = f(pair, output)
            end
        end
    end
    return output
end
