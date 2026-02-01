# ParticleSystem interface

The `ParticleSystem` interface facilitates the use of `CellListMap` for the majority of cases.

!!! note
    - This interface requires `CellListMap.jl` version `0.8.30` or greater.
    - The complete codes of the examples are available in the [Complete examples](@ref) section.

!!! compat
    The `ParticleSystem` interface is available since version `0.9.0` of CellListMap.jl.
    It replaces the `PeriodicSystems` interface available in previous versions.

## The mapped function

The purpose of CellListMap is to compute a pairwise-dependent function for all pairs of particles that are closer to
each other than a defined cutoff. This pairwise function must be implemented by the user and adhere to the following
interface:

```julia
function f(x, y, i, j, d2, output)
    # update output variable
    return output
end
```
where `x` and `y` are the positions of the particles, already wrapped relative to each other according to the periodic boundary conditions (a minimum-image set of positions), `i` and `j` are the indexes of the particles in the arrays of coordinates, `d2` is the squared distance between the particles, and `output` is the variable to be computed.

!!! info
    #### Details of the mapped function interface

    The input parameters `x`, `y`, `i`, `j`, and `d2` must not be modified by the user. They are the
    the input data that the user may use to update the `output` variable.

    |  Input Parameter |  Type   |  Meaning   |
    |:----------------:|:-------:|:-----------|
    | `x`              |  `SVector`   |  The coordinates of particle `i` of the pair.  |
    | `y`              |  `SVector`   |  The coordinates of particle `j` of the pair (minimum-image relative to `x`).  |
    | `i`              |  `Int`       |  Index of first particle in the original array of coordinates. |
    | `j`              |  `Int`       |  Index of second particle in the original array of coordinates. |
    | `d2`             |  `<:Real`    |  Squared distance between the particles. |
    | `output`         |  user defined   |  the value to be updated  |

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
function energy(d2,u)
    u += 1 / sqrt(d2)
    return u
end
```
and the additional parameters required by the interface can be eliminated by the use of an anonymous function, directly on the call to the `map_pairwise` function:
```julia
u = map_pairwise(
    (x,y,i,j,d2,u) -> energy(d2,u),
    system
)
```
(what `system` is will be explained in the examples below). Note that the `energy` function does not use the `x`, `y`, `i`, and `j` input parameters, such
that the anonymous function managing the interface could also be written as `(_, _, _, _, d2, u) -> energy(d2, u)`, making explicit the dummy character of
these variables in the example.

Alternatively, the function might require additional parameters, such as the masses of the particles. In this case, we can use a closure to provide such data:
```julia
function energy(i,j,d2,u,masses)
    u += masses[i]*masses[j] / sqrt(d2)
    return u
end
const masses = # ... some masses
u = map_pairwise((x,y,i,j,d2,u) -> energy(d2,u,masses), system)
```

Here we reinforce the fact that the `energy` functions defined above compute the contribution to the energy of the interaction of *a single* pair
of particles. This function will be called for every pair of particles within the cutoff, automatically, in the `map_pairwise` call.

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
