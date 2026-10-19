#
# Wrapper of the list of neighbors, that allows in-place updating of the lists
#
mutable struct NeighborList{T}
    n::Int
    list::Vector{Tuple{Int, Int, T}}
end

# resize! used to sizehint the neighborlist in a 1.10 compatible way.
# (no shrink=false)
function Base.resize!(nb::NeighborList, n) 
    resize!(nb.list, n)
    return nb
end

# Functions for parallel construction and reduction
function reset_output!(nb::NeighborList)
    nb.n = 0
    return nb
end
copy_output(nb::NeighborList) = NeighborList(nb.n, copy(nb.list))

# Specialized reduce to avoid O(nbatches²) work from incremental resize! calls:
# compute total once, resize output.list once, then copyto! from each batch.
# Large lists are copied in parallel (each batch to its own range of the output).
const PARALLEL_REDUCE_MIN_PAIRS = 2^16
function reduce_output!(output::NeighborList, output_threaded::Vector{<:NeighborList})
    ntot = output.n
    for nb in output_threaded
        ntot += nb.n
    end
    if length(output.list) < ntot
        resize!(output.list, ntot)
    end
    if nthreads() > 1 && length(output_threaded) > 1 && ntot - output.n >= PARALLEL_REDUCE_MIN_PAIRS
        tasks = Vector{Task}(undef, length(output_threaded))
        offset = output.n
        for (k, nb) in enumerate(output_threaded)
            tasks[k] = @spawn copyto!($output.list, $offset + 1, $nb.list, 1, $nb.n)
            offset += nb.n
        end
        foreach(wait, tasks)
    else
        offset = output.n
        for nb in output_threaded
            copyto!(output.list, offset + 1, nb.list, 1, nb.n)
            offset += nb.n
        end
    end
    output.n = ntot
    return output
end

# Estimate the expected number of neighbor pairs assuming uniform density.
function _estimated_n_pairs(box::Box{<:UnitCellType, N}, cl::CellList) where {UnitCellType, N}
    n = cl.n_real_particles
    vol = abs(LinearAlgebra.det(box.input_unit_cell.matrix))
    sphere_vol = (N == 2 ? π : 4π / 3) * box.cutoff^N
    return max(0, round(Int, (n * (n - 1) / 2) * sphere_vol / vol))
end
function _estimated_n_pairs(box::Box{<:UnitCellType, N}, cl::CellListPair) where {UnitCellType, N}
    n_x = cl.ref_list.n_real_particles
    n_y = cl.target_list.n_real_particles
    vol = abs(LinearAlgebra.det(box.input_unit_cell.matrix))
    sphere_vol = (N == 2 ? π : 4π / 3) * box.cutoff^N
    return max(0, round(Int, n_x * n_y * sphere_vol / vol))
end

# Pre-allocate capacity for neighbor lists based on the estimated pair count: the
# output list holds all pairs, and each threaded list holds the pairs of one batch
# (with a margin for load imbalance). Resizing to a smaller length keeps the
# capacity, such that repeated calls do not allocate.
function _sizehint_neighbor_lists!(sys::AbstractParticleSystem)
    n_pairs = _estimated_n_pairs(box(sys), celllist(sys))
    resize!(sys.output, n_pairs)
    _output_threaded = output_threaded(sys)
    n_per_batch = min(n_pairs, cld(n_pairs, max(1, length(_output_threaded))) * 5 ÷ 4)
    for nb in _output_threaded
        resize!(nb, max(length(nb.list), n_per_batch))
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

