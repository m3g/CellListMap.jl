# Migrating from 0.9

This guide helps users migrate from CellListMap version 0.9.x to version 0.10.0.

## Neighbor lists

The `neighborlist` and `neighborlist!` functions are now called with keyword parameters following the same interface of `ParticleSystem`. This means:

```julia
# Before (0.9.x)
neighborlist(x, 0.01; unitcell=[1,1,1])

# After (0.10.x)
neighborlist(positions=x, cutoff=0.01, unitcell=[1,1,1])

# Before (0.9.x)
neighborlist(x, y, 0.01; unitcell=[1,1,1])

# After (0.10.x)
neighborlist(xpositions=x, ypositions=y, cutoff=0.01, unitcell=[1,1,1])
```

## InPlaceNeighborList

Similarly, `InPlaceNeighborList` constructor follows the same new keyword syntax:

```julia
# Before (0.9.x)
InPlaceNeighborList(x, 0.01; unitcell=[1,1,1])

# After (0.10.x)
InPlaceNeighborList(positions=x, cutoff=0.01, unitcell=[1,1,1])

# Before (0.9.x)
InPlaceNeighborList(x, y, 0.01; unitcell=[1,1,1])

# After (0.10.x)
InPlaceNeighborList(xpositions=x, ypositions=y, cutoff=0.01, unitcell=[1,1,1])
```

The `update!` function was similarly updated, with the same keyword arguments.

### Removal of `autoswap`

The `autoswap` option, which was deprecated in version 0.9.16, has been completely removed in 0.10.0. This option was previously used to automatically swap the sets of particles to optimize performance when computing neighbor lists between two sets of particles with different sizes.

**Migration**: If you were using `autoswap=true` (which was the default), you may experience some performance regression for smaller systems. For most use cases, simply remove the `autoswap` keyword argument.

## ParticleSystem interface updates

### Renaming of `map_pairwise!` to `pairwise!`

The function `map_pairwise!` has been renamed to `pairwise!`:

```julia
# Before (0.9.x)
u = map_pairwise!(f, system)

# After (0.10.0)
u = pairwise!(f, system)
```

### Removal of `pairwise` (without `!`)

The non-mutating `pairwise` function has been removed. This change emphasizes that the function always mutates the `output` field of the `ParticleSystem` object:

```julia
# Before (0.9.x)
u = pairwise(f, system)

# After (0.10.0)
u = pairwise!(f, system)  # always use the mutating version
```

### Simplified mapped function signature

The signature of the mapped function has been simplified. Instead of receiving individual arguments, the function now receives a `NeighborPair` object containing all pair information:

```julia
# Before (0.9.x)
function f(x, y, i, j, d2, output)
    d = sqrt(d2)
    # ... compute something
    return output
end

# After (0.10.0)
function f(pair, output)
    (; x, y, i, j, d2, d) = pair  # destructuring
    # ... compute something
    return output
end
```

The `NeighborPair` object contains the fields:
- `x`, `y`: positions of the two particles (minimum-image convention applied)
- `i`, `j`: indices of the particles in the original arrays
- `d2`: squared distance between the particles
- `d`: distance between the particles (**lazily computed** - only calculated when accessed)

The lazy computation of `d` is particularly useful when your function only needs the squared distance `d2`, avoiding unnecessary `sqrt` calls.

### Output resetting behavior

The handling of the initial output value has changed:

- **Before (0.9.x)**: The `output` was reset to zero always, and the initial value given to `ParticleSystem` was ignored. 
- **After (0.10.0)**: The `output` is reset to `zero(typeof(output))` by default at each `pairwise!` call, but the `reset=false` option can be used to accumulate with the initial value.

To skip the automatic resetting (useful for accumulating results), use the `reset=false` keyword:

```julia
# Reset output to zero before computation (default)
pairwise!(f, system)

# Keep current output value and accumulate
pairwise!(f, system; reset=false)
```

## `ParticleSystem` updating

The interface for updating particle system was simplified. `update_cutoff!`, `update_unicell!` and direct access to the positions was removed. Now, all updates must be done through the `update!` function as, for example:

```julia
update(system; 
    cutoff=0.02, p
    positions=rand(3,1000), # or xpositions=
    ypositions=rand(3,1000), # if ypositions were defined
    unitcell=[1,1,1,2], 
    parallel=false
)
```

## Low-level interface

The previous "low-level interface" using `Box`, `CellList`, and `pairwise!(f, output, box, cl)` directly is now internal and no longer exported or part of the public API.

**Migration**: Use the `ParticleSystem` interface instead. The `ParticleSystem` interface provides all the functionality of the low-level interface with a simpler and more user-friendly API.

```julia
# Before (0.9.x) - Low-level interface
box = Box(sides, cutoff)
cl = CellList(x, box)
u = map_pairwise!((pair, u) -> f(pair, u), 0.0, box, cl)

# After (0.10.0) - ParticleSystem interface
sys = ParticleSystem(positions=x, cutoff=cutoff, unitcell=sides, output=0.0)
u = pairwise!(f, sys)
```

For iterative computations where the particle positions change:

```julia
# Before (0.9.x)
box = Box(sides, cutoff)
cl = CellList(x, box)
for step in 1:nsteps
    x = new_positions()
    UpdateCellList!(x, box, cl)
    u = map_pairwise!(f, 0.0, box, cl)
end

# After (0.10.0)
sys = ParticleSystem(positions=x, cutoff=cutoff, unitcell=sides, output=0.0)
for step in 1:nsteps
    update!(sys; xpositions=new_positions())
    u = pairwise!(f, sys)
end
```

See the [ParticleSystem](@ref ParticleSystem-interface) and [Updating the system](@ref) sections for more details on using the `ParticleSystem` interface.

## Python interface

The Python interface has been discontinued in version 0.10.0.
