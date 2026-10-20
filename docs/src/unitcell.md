# Unit cell size and cutoff

There is no restriction on the size of the unit cell relative to the cutoff: the cutoff
may even be larger than the unit cell itself. The only requirement is that the unit cell
is not degenerate (it must have a non-zero volume), otherwise an error is thrown.

The mapped function of `pairwise!` is called once for each periodic image of particle `j`
that is within the cutoff of particle `i`, and neighbor lists contain one entry for each
such image. In the most common situation, the unit cell is larger than twice the cutoff
(for orthorhombic cells, each side is greater than twice the cutoff; for triclinic cells,
the distances between opposite faces are). Then, at most one image of each particle can be within
the cutoff of another particle, each pair of particles is reported at most once, and
`pair.y` is the minimum image of particle `j` relative to `pair.x`.

## Small unit cells

If the unit cell is smaller (for example, in crystals with small unit cells, or
when using long cutoffs), more than one periodic image of a particle may be within the cutoff
of another particle, and a particle may be within the cutoff of its own images. Thus:

- the same pair `(i, j)` may appear more than once, with different distances; 
- pairs with `i == j` may appear, corresponding to a particle and one of its own images.

This is the physically correct behavior for the computation of properties of periodic systems,
as for example the energy of a crystal: the contributions of all images within the cutoff are summed.

For example, a single particle in a cubic cell of side `1.0`, with a cutoff of `1.5`, is
within the cutoff of 9 of its own images (the 6 face-neighbors at distance 1 and the 12
edge-neighbors at distance √2, each pair being reported once):

```@example unitcell
using CellListMap, StaticArrays
x = [SVector(0.5, 0.5, 0.5)]
neighborlist(; positions = x, cutoff = 1.5, unitcell = [1.0, 1.0, 1.0])
```

!!! note
    The cost of the computation increases with the number of images within the cutoff. For cells much smaller
    than the cutoff, the number of images of each particle within the cutoff grows with $(r_{cut}/L)^3$, where $L$ is the size of the cell.

## Identifying the image: `pair.shift`

The `shift` field of the `NeighborPair` object passed to the mapped function identifies the image
of particle `j` in each pair: it is the integer lattice translation of particle `j` relative to particle `i`, such that

```julia
y[pair.j] + unitcell * pair.shift - x[pair.i] ≈ pair.y - pair.x
```

where `x` and `y` are the input coordinates (with `y = x` for single-set computations), and `unitcell` is the
unit cell matrix (or `Diagonal(sides)` for orthorhombic cells). The shift is defined relative to the input coordinates,
thus it also accounts for the wrapping of coordinates that are outside the unit cell. For
non-periodic systems, it is always zero.
