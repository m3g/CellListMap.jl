# Single set: Simple outputs

This section shows examples of computing simple outputs (a single scalar or a single array)
for a single set of particles.

## Potential energy example

For example, here we read the coordinates of Argon atoms from a PDB file. The coordinates are given as
vector of `SVector`s. We then compute an "energy", which in this case is simply the sum of `1/d` over all pair of particles, within a cutoff.

```@example energy
using CellListMap, PDBTools
argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))
system = ParticleSystem(
    xpositions=argon_coordinates,
    unitcell=[21.0,21.0,21.0],
    cutoff = 8.0,
    output = 0.0,
    output_name = :energy
)
```

Now, let us compute the energy of the particles, assuming a simple formula which depends on the inverse of the distance between pairs:

```@example energy
function energy(pair, energy)
    energy += 1 / pair.d
    return energy
end
```

```@example energy
pairwise!(energy, system)
```

Because `output_name` was set to `:energy`, the `system.energy` field accesses the resulting value of the computation:
```@example energy
system.energy
```
If the `output_name` field is not provided, the output value from the `system.output` field.

## The initial value of output

The `output = 0.0` in the construction of the example above is used, by default,
to set the type of the output variable, here a `Float64`. This value is stored in 
the `output` field, but it is, by default, reset to `zero(typeof(output))` when
calling the `pairwise!` function.

To use the initial value provided and *accumulate* on top of it in the call to `pairwise!`, 
the `reset=false` option must be provided, as:

```@example energy
pairwise!(energy, system; reset=false)
```

Note that the result is twice the previous value, because the initial value of the energy was not reset.
If we instead call `pairwise!` with the default parameters, we recompute the energy from scratch:

```@example energy
pairwise!(energy, system) # reset=true by default
```

## Computing forces

Following the example above, let us compute the forces between the particles. We have to define the function that computes the force between a pair of particles and updates the array of forces:

```@example energy
function update_forces!(pair, forces)
    (; i, j, x, y, d2, d) = pair
    df = (1/d2)*(1/d)*(y - x)
    forces[i] += df
    forces[j] -= df
    return forces
end
```

Importantly, the function *must* return the `forces` array to follow the API.

Now, let us setup the system with the new type of output variable, which will be now an array of forces with the same type as the positions:

```@example energy
using PDBTools
argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))
system = ParticleSystem(
    xpositions=argon_coordinates,
    unitcell=[21.0, 21.0, 21.0],
    cutoff = 8.0,
    output = similar(argon_coordinates),
    output_name = :forces
)
```

A call to `pairwise!` with the appropriate function definition will update the forces:
```@example energy
pairwise!((pair, forces) -> update_forces!(pair, forces), system)
```

If we want now to accumulate additional forces, we need to use `reset=false` in the call to `pairwise!`. 
