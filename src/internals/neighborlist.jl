#
# Wrapper of the list of neighbors, that allows in-place updating of the lists
#
mutable struct NeighborList{T}
    n::Int
    list::Vector{Tuple{Int, Int, T}}
end

# Functions for parallel construction and reduction
function reset_output!(nb::NeighborList)
    nb.n = 0
    return nb
end
copy_output(nb::NeighborList) = NeighborList(nb.n, copy(nb.list))
function reducer!(nb1::NeighborList, nb2::NeighborList)
    ntot = nb1.n + nb2.n
    if length(nb1.list) < ntot
        resize!(nb1.list, ntot)
    end
    nb1.list[nb1.n + 1:ntot] .= @view(nb2.list[1:nb2.n])
    nb1.n = ntot
    return nb1
end

# Specialized reduce to avoid O(nbatches²) work from incremental resize! calls:
# compute total once, resize output.list once, then copyto! from each batch.
function reduce_output!(output::NeighborList, output_threaded::Vector{<:NeighborList})
    ntot = output.n
    for nb in output_threaded
        ntot += nb.n
    end
    if length(output.list) < ntot
        resize!(output.list, ntot)
    end
    offset = output.n
    for nb in output_threaded
        output.list[offset + 1:offset + nb.n] .= @view(nb.list[1:nb.n])
        offset += nb.n
    end
    output.n = ntot
    return output
end

# Estimate the expected number of neighbor pairs assuming uniform density.
function _estimated_n_pairs(box::Box{UnitCellType, N}, cl::CellList) where {UnitCellType, N}
    n = cl.n_real_particles
    vol = abs(LinearAlgebra.det(box.input_unit_cell.matrix))
    sphere_vol = (N == 2 ? π : 4π / 3) * box.cutoff^N
    return max(0, round(Int, (n * (n - 1) / 2) * sphere_vol / vol))
end
function _estimated_n_pairs(box::Box{UnitCellType, N}, cl::CellListPair) where {UnitCellType, N}
    n_x = cl.ref_list.n_real_particles
    n_y = cl.target_list.n_real_particles
    vol = abs(LinearAlgebra.det(box.input_unit_cell.matrix))
    sphere_vol = (N == 2 ? π : 4π / 3) * box.cutoff^N
    return max(0, round(Int, n_x * n_y * sphere_vol / vol))
end

# Pre-allocate capacity for neighbor lists based on the estimated pair count.
# shrink=false ensures this is a no-op when capacity is already sufficient,
# preventing spurious allocations on repeated calls.
function _sizehint_neighbor_lists!(output::NeighborList, output_threaded, box, cl)
    n_pairs = _estimated_n_pairs(box, cl)
    sizehint!(output.list, n_pairs; shrink = false)
    nbatch = max(1, nbatches(cl, :map))
    n_per_batch = cld(n_pairs, nbatch)
    for nb in output_threaded
        sizehint!(nb.list, n_per_batch; shrink = false)
    end
    return nothing
end

# Function adds pair to the list
function push_pair!(pair, nb::NeighborList)
    (; i, j, d) =  pair
    nb.n += 1
    if nb.n > length(nb.list)
        push!(nb.list, (i, j, d))
    else
        nb.list[nb.n] = (i, j, d)
    end
    return nb
end

