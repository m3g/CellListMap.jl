# Single set: Simple outputs

This section shows examples of computing simple outputs (a single scalar or a single array)
for a single set of particles.

## Potential energy example

For example, here we read the coordinates of Argon atoms from a PDB file. The coordinates are given as
vector of `SVector`s. We then compute an "energy", which in this case is simply the sum of `1/d` over all pair of particles, within a cutoff.

```julia-repl
julia> using CellListMap, PDBTools

julia> argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))

julia> system = ParticleSystem(
           xpositions=argon_coordinates,
           unitcell=[21.0,21.0,21.0],
           cutoff = 8.0,
           output = 0.0,
           output_name = :energy
       );
```

Now, let us compute the energy of the particles, assuming a simple formula which depends on the inverse of the distance between pairs:

```julia-repl
function energy(pair, energy)
    energy += 1 / pair.d
    return energy
end

julia> map_pairwise!(energy, system)
207.37593043370865
```
Note that the first four parameters of `energy` are not used here but are needed to adhere to the interface. The function
input could be written as `(_, _, _, _, d2, energy)` to make that explicit.

Because `output_name` was set to `:energy`, the `system.energy` field accesses the resulting value of the computation:
```julia-repl
julia> system.energy
207.37593043370865
```
If the `output_name` field is not provided, the output value from the `system.output` field.

## Computing forces

Following the example above, let us compute the forces between the particles. We have to define the function that computes the force between a pair of particles and updates the array of forces:

```julia
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

```julia-repl
julia> using CellListMap, PDBTools

julia> argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))

julia> system = ParticleSystem(
           xpositions=argon_coordinates,
           unitcell=[21.0, 21.0, 21.0],
           cutoff = 8.0,
           output = similar(argon_coordinates),
           output_name = :forces
       );
```

Let us note that the `forces` where reset upon the construction of the system:
```julia-repl
julia> system.forces
1000-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 0.0]
 ⋮
 [0.0, 0.0, 0.0]
```

A call to `map_pairwise!` with the appropriate function definition will update the forces:
```julia-repl
map_pairwise!((pair, forces) -> update_forces!(pair, forces), system)
100-element Vector{SVector{3, Float64}}:
 [0.026493833307357332, 0.18454277989323772, -0.012253902366284965]
 [0.07782602581235695, 0.2791082233740261, 0.21926615329195248]
 ⋮
 [0.11307234751448932, 0.006353545239676281, -0.05955687310348302]
 [-0.03101200918307673, 0.03543655648545697, 0.031849121630976335]
```
