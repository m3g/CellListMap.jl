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
if the number of particles changes.

### Example: updating in a simulation loop

```julia-repl
julia> using CellListMap, StaticArrays

julia> system = ParticleSystem(
           xpositions = rand(SVector{3,Float64}, 1000),
           unitcell = [1.0, 1.0, 1.0],
           cutoff = 0.1,
           output = 0.0,
       );

julia> new_x = rand(SVector{3,Float64}, 1000);

julia> update!(system; xpositions=new_x, unitcell=SVector(1.2, 1.2, 1.2), cutoff=0.15);
```

Multiple properties can be updated in a single call. Any combination of kwargs is valid.

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

### Fine-grained mutation via `system.xpositions`

Individual positions can still be updated directly through the `xpositions` (and
`ypositions`) property, which behaves as an ordinary `Vector{SVector{N,T}}`:

```julia-repl
julia> system.xpositions[1] = SVector(0.0, 0.0, 0.0)

julia> push!(system.xpositions, SVector(0.5, 0.5, 0.5))

julia> system.xpositions .= new_positions   # in-place broadcast
```

All of these operations set an internal flag so that cell lists are recomputed
on the next call to `pairwise!`. The `output` variable may need to be resized
with `resize_output!` if the number of particles changed.

## Updating the unit cell

```julia-repl
julia> update!(system; unitcell=SVector(1.2, 1.2, 1.2))                          # orthorhombic

julia> update!(system; unitcell=[1.2 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.2])        # triclinic
```

The property setter `system.unitcell = x` is also available as a shortcut.

!!! note
    - The unit cell type (orthorhombic or triclinic) cannot be changed after construction.
    - The unit cell of non-periodic systems (initialized with `nothing`) cannot be updated.
    - Using static arrays or `Tuple`s makes the update non-allocating.

## Updating the cutoff

```julia-repl
julia> update!(system; cutoff=0.2)
```

The property setter `system.cutoff = 0.2` is also available as a shortcut.

## Updating the parallelization flag

```julia-repl
julia> update!(system; parallel=false)   # disable multi-threading

julia> update!(system; parallel=true)    # enable multi-threading
```
