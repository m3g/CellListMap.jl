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

```@docs
NeighborPair
ParticleSystem
```

### The parwise! methods

```@docs
pairwise!(::F, ::CellListMap.AbstractParticleSystem) where {F<:Function} 
pairwise!(::F, ::AbstractVecOrMat, ::CellListMap.ParticleSystem1) where {F<:Function}
```

### Updating systems

```@docs
update_cutoff!
update_unitcell!
resize_output!
```

### Custom parallel reduction

These are public, but not exported.

```@docs
CellListMap.copy_output
CellListMap.reset_output!
CellListMap.reducer!
```

### Auxiliary functions

These are public, but not exported.

```@docs
CellListMap.wrap_relative_to
CellListMap.get_computing_box
```
