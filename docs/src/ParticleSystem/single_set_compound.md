# Single set: Compound outputs

This section shows how to compute compound outputs (e.g., both energy and forces)
in a single computation pass, which requires defining custom output types.

## Computing both energy and forces

In this example we define a general type of `output` variable, for which custom copy, reset, and reduction functions
must be defined. It can be followed for the computation of other general properties from the particle positions.

!!! note
    Interface to be implemented:

    |   Method              |   Return    |  What it does    |
    |:----------------------|:-----------:|:---------------- |
    | `copy_output(x::T)`   | new instance of type `T` | Copies an element of the output type `T`. |
    | `reset_output!(x::T)` | mutated `x` | Resets (usually zero) the value of x to the initial value it must assume before mapping.  If `x` is immutable, the function can return a new instance of `T`. |
    | `reducer(x::T,y::T)` | `mutated x` | Reduces `x` and `y` into `x` (for example `x = x + y`). If `x` is immutable, returns a new instance of type `T`.

    **Remark:** if the output is an array of an immutable type `T`, the methods above can be defined for single *instances* of `T`, which is simpler
    than for the arrays.

```julia
using CellListMap, StaticArrays, PDBTools
```

The computation of energies and forces in a single call is an interesting example for the definition of a custom `output` type and the required interface functions.
Let us first define an output variable containing both quantities:
```julia
mutable struct EnergyAndForces
    energy::Float64
    forces::Vector{SVector{3,Float64}}
end
```

Now we need to define what it means to copy, reset, and reduce this new type of output. We overload
the default corresponding functions, for our new output type:

The copy method creates a new instance of the `EnergyAndForces` type, with copied data:
```julia
function CellListMap.copy_output(x::EnergyAndForces)
    return EnergyAndForces(copy(x.energy), copy(x.forces))
end
```

The reset method will zero both the energy and all forces:
```julia
function CellListMap.reset_output!(output::EnergyAndForces)
    output.energy = 0.0
    for i in eachindex(output.forces)
        output.forces[i] = SVector(0.0, 0.0, 0.0)
    end
    return output
end
```

The reducer function defines what it means to combine two output variables obtained on
independent threads. In this case, we sum the energies and forces. Different reduction functions
might be necessary for other custom types (for example if computing minimum distances).
```julia
function CellListMap.reducer(x::EnergyAndForces, y::EnergyAndForces)
    e_tot = x.energy + y.energy
    x.forces .+= y.forces
    return EnergyAndForces(e_tot, x.forces)
end
```
Note that in the above example, we reuse the `x.forces` array in the return instance
of `EnergyAndForces`. You must always reduce from right to left, and reuse the
possible buffers of the first argument of the reducer (in this case, `x`).

!!! warning
    - All these functions **must** return the modified `output` variable, to adhere to the interface.
    - The proper definition of a reduction function is crucial for correctness. Please verify
      your results if using the default reducer function, which sums the elements.

Now we can proceed as before, defining a function that updates the output variable appropriately:
```julia
function energy_and_forces!(pair, output::EnergyAndForces)
    d = pair.d
    output.energy += 1/d
    df = (1/pair.d2)*(1/d)*(pair.y - pair.x)
    output.forces[pair.i] += df
    output.forces[pair.j] -= df
    return output
end
```

To finally define the system and compute the properties:

```julia
argon_coordinates = coor(read_pdb(CellListMap.argon_pdb_file))

system = ParticleSystem(
    xpositions = argon_coordinates,
    unitcell = [21.0, 21.0, 21.0],
    cutoff = 8.0,
    output = EnergyAndForces(0.0, similar(argon_coordinates)),
    output_name = :energy_and_forces
);

foreachneighbor!(energy_and_forces!, system);
```

The output can be seen with the aliases of the `system.output` variable:
```julia-repl
julia> system.energy_and_forces.energy
207.37593043370862

julia> system.energy_and_forces.forces
100-element Vector{SVector{3, Float64}}:
 [0.02649383330735732, 0.18454277989323772, -0.012253902366284958]
 [0.07782602581235692, 0.27910822337402613, 0.21926615329195248]
 ⋮
 [0.11307234751448932, 0.006353545239676281, -0.05955687310348303]
 [-0.031012009183076745, 0.03543655648545698, 0.03184912163097636]
```
