#
# Pairwise computations: loops over the pairs of particles within the cutoff,
# for self (`CellList`) and cross (`CellListPair`) computations.
#
# For each cell containing particles, the neighboring cells are visited through
# the precomputed offsets. For periodic systems, the index of the neighboring cell
# is wrapped (`divmod`), and the quotient is the lattice translation of that
# cell, which brings the image of its particles close to the particles of the
# current cell. For non-periodic systems, offsets that fall outside the grid
# are skipped.
#

#=
    _pairwise!(
        f::Function,
        output,
        box::Box,
        cl::Union{CellList,CellListPair};
        parallel::Bool=true,
        output_threaded=nothing,
        show_progress::Bool=false
    )

This function will run over every pair of particles which are closer than
`box.cutoff` and compute the Euclidean distance between the particles,
considering the periodic boundary conditions given in the `Box` structure.
If the distance is smaller than the cutoff, a function `f` of the
coordinates of the two particles will be computed.

The function `f` receives two arguments as input:
    f(pair, output)
Where `pair` is a `NeighborPair` struct with fields:
- `pair.i`, `pair.j`: indices of the particles
- `pair.x`, `pair.y`: coordinates (of the image of `j` within the cutoff of `i`)
- `pair.d`: distance between particles
- `pair.d2`: squared distance
- `pair.shift`: lattice translation of the image of `j` relative to `i`

The function must return the updated `output` variable.

=#
function _pairwise!(
        f::F, output, box::Box, cl::Union{CellList, CellListPair};
        parallel::Bool = true,
        output_threaded = nothing,
        show_progress::Bool = false,
    ) where {F} # Needed for specialization for this function (avoids some allocations)
    if parallel
        output = _pairwise_parallel!(f, output, box, cl; output_threaded, show_progress)
    else
        output = _pairwise_serial!(f, output, box, cl; show_progress)
    end
    return output
end

# The list whose cells are iterated in the outer loop
_outer_list(cl::CellList) = cl
_outer_list(cl::CellListPair) = cl.ref_list

function _pairwise_serial!(
        f::F, output, box::Box, cl::Union{CellList, CellListPair};
        show_progress::Bool = false
    ) where {F}
    cells = _outer_list(cl).nonempty_cells
    p = show_progress ? Progress(length(cells), dt = 1) : nothing
    for c in cells
        output = _cell_pairs!(f, output, box, cl, c)
        _next!(p)
    end
    return output
end

#
# Work-queue scheduling: the range of cells is split into many small consecutive chunks,
# and a fixed pool of `_nbatches` worker tasks (one per `output_threaded` accumulator
# slot) repeatedly claims the next unclaimed chunk from a shared atomic counter, such that
# faster workers pull more chunks.
#
function _batch!(f::F, ibatch, chunks, next_chunk, output_threaded, box, cl, p) where {F}
    cells = _outer_list(cl).nonempty_cells
    while true
        ic = atomic_add!(next_chunk, 1)
        ic > length(chunks) && break
        for k in chunks[ic]
            output_threaded[ibatch] = _cell_pairs!(f, output_threaded[ibatch], box, cl, cells[k])
            _next!(p)
        end
    end
    return
end

function _pairwise_parallel!(
        f::F, output, box::Box, cl::Union{CellList, CellListPair};
        output_threaded = nothing,
        show_progress::Bool = false
    ) where {F}
    _nbatches = nbatches(cl, :map)
    if isnothing(output_threaded)
        output_threaded = [deepcopy(output) for _ in 1:_nbatches]
    end
    cells = _outer_list(cl).nonempty_cells
    p = show_progress ? Progress(length(cells), dt = 1) : nothing
    n_chunks = _n_workqueue_chunks(_nbatches, _outer_list(cl).n_real_particles)
    chunks = collect(index_chunks(1:length(cells); n = n_chunks, split = Consecutive()))
    next_chunk = Atomic{Int}(1)
    tasks = Vector{Task}(undef, _nbatches)
    for ibatch in 1:_nbatches
        tasks[ibatch] = @spawn _batch!($f, $ibatch, $chunks, $next_chunk, $output_threaded, $box, $cl, $p)
    end
    foreach(wait, tasks)
    return reduce_output!(output, output_threaded)
end

#
# Neighboring cell of the cell with (0-based) cartesian indices `ci`, through offset `o`.
# Returns a flag indicating if the cell exists, its linear index, and the lattice translation
# (in units of unit cells) of the neighboring cell.
#
@inline function _neighbor_cell(box::Box{<:PeriodicCellType, N}, ci::SVector{N, Int}, o::SVector{N, Int}) where {N}
    (; nc) = box
    t = ci + o
    q = fld.(t, nc)
    return true, _cell_linear_index(t - q .* nc, nc), q
end
@inline function _neighbor_cell(box::Box{NonPeriodicCell, N}, ci::SVector{N, Int}, o::SVector{N, Int}) where {N}
    (; nc) = box
    t = ci + o
    inside = all(t .>= 0) & all(t .< nc)
    return inside, (inside ? _cell_linear_index(t, nc) : 0), zero(SVector{N, Int})
end

@inline _shift(s::SVector{N, Int32}) where {N} = SVector{N, Int}(s)

#
# Pairs involving the particles of cell `c`, for self-computations: pairs within the cell
# (each pair once), and pairs with the particles of the neighboring cells reached by the
# half set of offsets (each pair of cells once, including the images of the cell itself
# in small unit cells).
#
function _cell_pairs!(f::F, output, box::Box{UnitCellType, N, T}, cl::CellList{N, T}, c::Int) where {F, UnitCellType, N, T}
    (; cutoff_sqr, nc) = box
    m = box.input_unit_cell.matrix
    particles = cl.particles
    shifts = cl.shifts
    cell_start = cl.cell_start
    @inbounds first_i, last_i = cell_start[c], cell_start[c + 1] - 1
    # Pairs within the same cell
    @inbounds for a in first_i:last_i
        pa = particles[a]
        xa = pa.coordinates
        for b in (a + 1):last_i
            pb = particles[b]
            xb = pb.coordinates
            d2 = sum(abs2, xb - xa)
            if d2 <= cutoff_sqr
                shift = _shift(shifts[b]) - _shift(shifts[a])
                output = f(NeighborPair(pa.index, pb.index, xa, xb, d2, shift), output)
            end
        end
    end
    # Pairs with the particles of neighboring cells
    ci = _cell_cartesian_index(c, nc)
    @inbounds for o in cl.offsets_half
        exists, cj, q = _neighbor_cell(box, ci, o)
        exists || continue
        first_j, last_j = cell_start[cj], cell_start[cj + 1] - 1
        first_j > last_j && continue
        cell_shift = m * q
        for a in first_i:last_i
            pa = particles[a]
            xa = pa.coordinates
            # shifting `xa` instead of `xb` such that the translation is applied once per particle
            xa_shifted = xa - cell_shift
            for b in first_j:last_j
                pb = particles[b]
                d2 = sum(abs2, pb.coordinates - xa_shifted)
                if d2 <= cutoff_sqr
                    shift = q + _shift(shifts[b]) - _shift(shifts[a])
                    output = f(NeighborPair(pa.index, pb.index, xa, pb.coordinates + cell_shift, d2, shift), output)
                end
            end
        end
    end
    return output
end

#
# Pairs involving the particles of cell `c` of the reference list, and the particles
# of the target list in the same and neighboring cells (all offsets).
#
function _cell_pairs!(f::F, output, box::Box{UnitCellType, N, T}, cl::CellListPair{N, T}, c::Int) where {F, UnitCellType, N, T}
    (; cutoff_sqr, nc) = box
    m = box.input_unit_cell.matrix
    ref = cl.ref_list
    target = cl.target_list
    @inbounds first_i, last_i = ref.cell_start[c], ref.cell_start[c + 1] - 1
    ci = _cell_cartesian_index(c, nc)
    @inbounds for o in ref.offsets_full
        exists, cj, q = _neighbor_cell(box, ci, o)
        exists || continue
        first_j, last_j = target.cell_start[cj], target.cell_start[cj + 1] - 1
        first_j > last_j && continue
        cell_shift = m * q
        for a in first_i:last_i
            pa = ref.particles[a]
            xa = pa.coordinates
            xa_shifted = xa - cell_shift
            for b in first_j:last_j
                pb = target.particles[b]
                d2 = sum(abs2, pb.coordinates - xa_shifted)
                if d2 <= cutoff_sqr
                    shift = q + _shift(target.shifts[b]) - _shift(ref.shifts[a])
                    output = f(NeighborPair(pa.index, pb.index, xa, pb.coordinates + cell_shift, d2, shift), output)
                end
            end
        end
    end
    return output
end
