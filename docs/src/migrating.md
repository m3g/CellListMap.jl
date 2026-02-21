# Migrating from 0.9

This guide helps users migrate from CellListMap version 0.9.x to version 0.10.0.

## Neighbor lists updates

The `neighborlist` and `neighborlist!` functions remain mostly unchanged in their interface, except for the removal of the `autoswap` option.

### Removal of `autoswap`

The `autoswap` option, which was deprecated in version 0.9.16, has been completely removed in 0.10.0. This option was previously used to automatically swap the sets of particles to optimize performance when computing neighbor lists between two sets of particles with different sizes.

**Migration**: If you were using `autoswap=true` (which was the default), you may experience some performance regression for smaller systems. For most use cases, simply remove the `autoswap` keyword argument:

```julia
# Before (0.9.x)
neighborlist(x, y, cutoff; autoswap=true)

# After (0.10.0)
neighborlist(x, y, cutoff)
```

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
    sys.xpositions .= new_positions()  # update coordinates directly
    u = pairwise!(f, sys)
end
```

See the [ParticleSystem](@ref ParticleSystem-interface) and [Updating the system](@ref) sections for more details on using the `ParticleSystem` interface.

## Python interface

The Python interface has been discontinued in version 0.10.0.

## Summary of breaking changes

| Change | 0.9.x | 0.10.0 |
|:-------|:------|:-------|
| Pairwise mapping function | `map_pairwise!` | `pairwise!` |
| Alias for "non-mutating" pairwise | `pairwise` | Removed |
| Output reset behavior | Manual | Automatic (use `reset=false` to skip) |
| Low-level interface | Exported | Internal (use `ParticleSystem`) |
| `autoswap` option | Deprecated | Removed |
| Python interface | Available | Discontinued |
