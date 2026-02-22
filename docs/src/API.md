```@meta
CollapsedDocStrings = true
```

# Application interface 

## Neighborlist interface

### Simple neighborlists

```@docs
neighborlist(::Any, ::Any)
neighborlist(::Any, ::Any, ::Any)
```

### In-place neighborlists

```@docs
InPlaceNeighborList
update!
neighborlist!
```

## ParticleSystems

### Structures

Public and exported:

```@docs
NeighborPair
ParticleSystem
```

Public but not exported:

```@docs
CellListMap.AbstractParticleSystem
CellListMap.ParticleSystem1
CellListMap.ParticleSystem2
```

### The pairwise! methods

```@docs
pairwise!(::F, ::CellListMap.AbstractParticleSystem) where {F<:Function}
```

### Updating systems

```@docs
update!(::CellListMap.AbstractParticleSystem)
resize_output!
```

!!! note "Deprecated update functions"
    `update_cutoff!` and `update_unitcell!` are deprecated. Use `update!(sys; cutoff=...)` and
    `update!(sys; unitcell=...)` instead. Setting `sys.parallel` directly is also deprecated;
    use `update!(sys; parallel=...)`.

```@docs
update_cutoff!
update_unitcell!
```

### Custom parallel reduction

These are public, but not exported.

```@docs
CellListMap.copy_output
CellListMap.reset_output!
CellListMap.reducer!
CellListMap.reduce_output!
```

### Auxiliary functions

These are public, but not exported.

```@docs
CellListMap.wrap_relative_to
CellListMap.get_computing_box
```
