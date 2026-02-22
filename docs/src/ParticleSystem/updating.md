# Updating the system

If the `pairwise!` function will compute energy and/or forces in an iterative procedure (a simulation, for instance), we need to update the coordinates, and perhaps the unit cell, the cutoff, or the parallelization flag.

## The `update!` function

All system properties can be updated with the single `update!` function. Only the keyword arguments that are provided are updated; the rest stay unchanged:

```julia
update!(system;
    xpositions = ...,   # new coordinates for the first set
    ypositions = ...,   # new coordinates for the second set (two-set systems only)
    cutoff     = ...,   # new cutoff distance
    unitcell   = ...,   # new unit cell
    parallel   = ...,   # true or false
)
```

Coordinates accept the same types as the `ParticleSystem` constructor: `Vector{SVector}`,
vector of plain vectors, or an `(D, N)` matrix. The internal storage is resized automatically
if the number of particles changes. When the number of particles does not change, coordinate
updates are non-allocating.

### Example: updating in a simulation loop

```julia
using CellListMap, StaticArrays
system = ParticleSystem(
    xpositions = rand(SVector{3,Float64}, 1000),
    unitcell = [1.0, 1.0, 1.0],
    cutoff = 0.1,
    output = 0.0,
)
for i in 1:nsteps
    new_x = rand(SVector{3,Float64}, 1000)
    update!(system; xpositions=new_x)
end
```

Multiple properties can be updated in a single call. Any combination of kwargs is valid.
Updating only coordinates, without changing the number of particles, does not incurr in new allocations. Updating the unitcell, the cutoff, or increasing the number of particles can incurr in new allocations because the structure of the cell lists can change. 

## Updating coordinates

Coordinates are updated by passing `xpositions` (and/or `ypositions` for two-set systems)
to `update!`. The input can be a `Vector{SVector}`, a vector of plain vectors, or a `(D,N)`
matrix — the same types accepted by the `ParticleSystem` constructor:

```julia-repl
julia> update!(system; xpositions=rand(SVector{3,Float64}, 1000))

julia> update!(system; xpositions=rand(3, 1000))              # (D,N) matrix

julia> update!(system; xpositions=[rand(3) for _ in 1:1000])  # vector-of-vectors
```

If the number of particles changes, the internal storage is resized automatically:

```julia-repl
julia> update!(system; xpositions=rand(SVector{3,Float64}, 1200))
# system now has 1200 particles
```

!!! warning
    The `output` variable may have to be resized accordingly when the number of particles
    changes. Use `resize_output!` (do **not** use `Base.resize!` directly):

    ```julia
    resize_output!(system, length(system.xpositions))
    ```

    In the case of compound outputs (custom output structures) like that
    of the [Computing both energy and forces](@ref) example, calling `resize_output!` will
    return an error and require the user to define `Base.resize!` for the
    custom type.

## Updating the unit cell

```julia
update!(system; unitcell=SVector(1.2, 1.2, 1.2)) # orthorhombic

update!(system; unitcell=[1.2 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.2]) # triclinic
```

!!! note
    - The unit cell type (orthorhombic or triclinic) cannot be changed after construction.
    - The unit cell of non-periodic systems (initialized with `nothing`) cannot be updated.

## Updating the cutoff

```julia-repl
julia> update!(system; cutoff=0.2)
```

## Updating the parallelization flag

```julia
update!(system; parallel=false) # disable multi-threading

update!(system; parallel=true) # enable multi-threading
```
