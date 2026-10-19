# Unit cell size and cutoff

There is no restriction on the size of the unit cell relative to the cutoff: the cutoff
may be larger than the unit cell itself. The only requirement is that the unit cell
is not degenerate (it must have a non-zero volume), otherwise an error is thrown.

## Large cells: one image per pair

In the most common situation, the distance between opposite faces of the unit cell is
greater than twice the cutoff. In this case, at most one periodic image of each particle
can be within the cutoff of any other particle (the minimum-image convention applies), and
each pair of particles is reported at most once, with `y` being the minimum image of
particle `j` relative to `x`.

For orthorhombic cells this means that each side is greater than twice the cutoff. In general, given
a unit cell matrix $\left[\vec{a}~\vec{b}~\vec{c}\right]$, where $\vec{a}$, $\vec{b}$, and $\vec{c}$ are
the lattice vectors, the conditions are:

```math
\vec{a} \cdot \frac{\vec{b} \times \vec{c}}{|\vec{b} \times \vec{c}|} \gt 2r_{cut}
```
```math
\vec{b} \cdot \frac{\vec{c} \times \vec{a}}{|\vec{c} \times \vec{a}|} \gt 2r_{cut}
```
```math
\vec{c} \cdot \frac{\vec{a} \times \vec{b}}{|\vec{a} \times \vec{b}|} \gt 2r_{cut}
```

where $r_{cut}$ is the cutoff. In 2D, given the unit cell matrix $\left[\vec{a}~\vec{b}\right]$, the conditions are:

```math
\sqrt{|\vec{a}|^2 - \left(\vec{a}\cdot \frac{\vec{b}}{|\vec{b}|}\right)^2} \gt 2r_{cut} 
```
```math
\sqrt{|\vec{b}|^2 - \left(\vec{b}\cdot \frac{\vec{a}}{|\vec{a}|}\right)^2} \gt 2r_{cut} 
```

## Small cells: one entry per image

If the unit cell is smaller (for example, in crystals with small unit cells, or
when using long cutoffs), more than one periodic image of a particle may be within the cutoff
of another particle, and a particle may be within the cutoff of its own images.
In this case, the mapped function of `pairwise!` is called once for each *image* within the cutoff,
and the neighbor lists contain one entry for each image. That is:

- the same pair `(i, j)` may appear more than once, with different distances; 
- pairs with `i == j` may appear, corresponding to a particle and one of its own images.

This is the physically correct behavior for the computation of properties of periodic systems,
as for example the energy of a crystal: the contributions of all images within the cutoff are summed.

The `shift` field of the `NeighborPair` object passed to the mapped function identifies the image:
it is the integer lattice translation of particle `j` relative to particle `i`, such that

```julia
y[pair.j] + unitcell * pair.shift - x[pair.i] ≈ pair.y - pair.x
```

where `x` and `y` are the input coordinates (with `y = x` for single-set computations), and `unitcell` is the
unit cell matrix (or `Diagonal(sides)` for orthorhombic cells). The shift is defined relative to the input coordinates,
thus it also accounts for the wrapping of coordinates that are outside the unit cell. For
non-periodic systems, it is always zero.

For example, a single particle in a cubic cell of side `1.0`, with a cutoff of `1.5`, is
within the cutoff of 9 of its own images (the 6 face-neighbors at distance 1 and the 12
edge-neighbors at distance √2, each pair being reported once):

```@example unitcell
using CellListMap, StaticArrays
x = [SVector(0.5, 0.5, 0.5)]
neighborlist(; positions = x, cutoff = 1.5, unitcell = [1.0, 1.0, 1.0])
```

!!! note
    The cost of the computation increases with the number of images within the cutoff, as each image is 
    explicitly replicated in the cell lists. For cells much smaller than the cutoff, the
    number of images per particle grows with $(r_{cut}/L)^3$, where $L$ is the size of the cell. 
