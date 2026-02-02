# ParticleSystem interface

## The mapped function

The purpose of CellListMap is to compute a pairwise-dependent function for all pairs of particles that are closer to
each other than a defined cutoff. This pairwise function must be implemented by the user and adhere to the following
interface:

```julia
function f(pair, output)
    (; i, j, x, y, d2, d) = pair
    # update output variable using above pair data
    return output
end
```

where `pair` is a `NeighborPair` object, containing fields `i`, `j`, `x`, `y`, `d2` and `d`, which here we
extract from the object using [destructuring syntax](https://docs.julialang.org/en/v1/manual/functions/#Property-destructuring). 

`x` and `y` are the positions of the particles, already wrapped relative to each other according to the periodic boundary conditions (a minimum-image set of positions), `i` and `j` are the indexes of the particles in the arrays of coordinates, `d2` is the squared distance between the particles, `d` is the distance, and `output` is the variable to be computed. 

!!! info
    #### Details of the mapped function interface

    The `pair` input of the function contains the  data that the user may use to update the `output` variable.
    The `NeighborPair` object contains fields `x`, `y`, `i`, `j`, `d2` and the lazily computed `d`, meaning:

    |  Input/Output |  Type   |  Meaning   |
    |:----------------:|:-------:|:-----------|
    | `x`              |  `SVector`   |  The coordinates of particle `i` of the pair.  |
    | `y`              |  `SVector`   |  The coordinates of particle `j` of the pair (minimum-image relative to `x`).  |
    | `i`              |  `Int`       |  Index of first particle in the original array of coordinates. |
    | `j`              |  `Int`       |  Index of second particle in the original array of coordinates. |
    | `d2`             |  `<:Real`    |  Squared distance between the particles. |
    | `d`             |  `<:Real`    |  Squared distance between the particles (computed lazily). |
    | `output`         |  user defined   |  the name of the variable to be updated.  |

    **Notes:** `x` and `y` may be 2D or 3D vectors, depending on the dimension of the system. The type of
    the coordinates of `x`, `y`, and of `d2` are dependent on the input arrays and cutoff, and can be `Float64`,
    `Float32`, unitful quantities, etc.

    |  Return value |  Type   |  Meaning   |
    |:----------------:|:-------:|:-----------|
    | `output`         |  user defined   |  the updated value of output.  |

    The `output` variable **must** be returned by the function, being it mutable or immutable.

### Basic examples

For example, computing the energy, as the sum of the inverse of the distance between particles, can be done with a function like:
```julia
function energy(pair, u)
    u += 1 / pair.d
    return u
end
```
and this function is passed directly to `pairwise!`:
```julia
u = pairwise!(energy, system)
```
(what `system` is will be explained in the examples below). Note that the `energy` function only uses `pair.d` (the distance), but all other fields (`pair.x`, `pair.y`, `pair.i`, `pair.j`, `pair.d2`) are available.

Alternatively, the function might require additional parameters, such as the masses of the particles. In this case, we can use a closure to provide such data:
```julia
function energy(pair, u, masses)
    (; i, j, d)  = pair 
    u += masses[i]*masses[j] / d
    return u
end
const masses = # ... some masses
u = pairwise!((pair, u) -> energy(pair, u, masses), system)
```

Here we reinforce the fact that the `energy` functions defined above compute the contribution to the energy of the interaction of *a single* pair
of particles. This function will be called for every pair of particles within the cutoff, automatically, in the `pairwise!` call.

!!! note
    The `output` of the `CellListMap` computation may be of any kind. Most commonly, it is an energy, a set of forces, or other data type that can be represented either as a number, an array of numbers, or an array of vectors (`SVectors` in particular), such as an arrays of forces.

    Additionally, the properties are frequently additive (the energy is the sum of the energy of the particles, or the forces are added by summation).

    For these types of `output` data the usage does not require the implementation of any data-type dependent function.

## The ParticleSystem constructor

The `ParticleSystem` constructor receives the properties of the system and sets up automatically the most commonly used data structures necessary.

!!! note
    - Systems can be 2 or 3-dimensional.
    - The `unitcell` parameter may be:
        - a vector, in which case the system periodic boundaries are Orthorhombic, this is faster.
        - a matrix, in which case the system periodic boundaries are Triclinic (general). The lattice
          vectors correspond to the *columns* of the matrix.
        - `nothing` (by default), in which case no periodic boundary conditions will be used.
    - `Unitful` quantities can be provided, given appropriate types for all input parameters.

## Docstrings

```@meta
CollapsedDocStrings = true
```

```@autodocs
Modules = [CellListMap]
Pages = ["ParticleSystem.jl"]
Order = [:function, :type]
```
