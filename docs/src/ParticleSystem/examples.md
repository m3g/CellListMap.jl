# Complete examples

This section contains complete, self-contained example codes that can be copied
and run directly.

## Simple energy computation

In this example, a simple potential energy defined as the sum of the
inverse of the distance between the particles is computed.

```julia
using CellListMap
using StaticArrays
system = ParticleSystem(
    xpositions = rand(SVector{3,Float64},1000),
    unitcell=[1.0,1.0,1.0],
    cutoff = 0.1,
    output = 0.0,
    output_name = :energy
)
map_pairwise!((pair, energy) -> energy += 1 / pair.d, system)
```

## Force computation

Here we compute the force vector associated to the potential energy
function of the previous example.

```julia
using CellListMap
using StaticArrays
positions = rand(SVector{3,Float64},1000)
system = ParticleSystem(
    xpositions = positions,
    unitcell=[1.0,1.0,1.0],
    cutoff = 0.1,
    output = similar(positions),
    output_name = :forces
)
function update_forces!(pair, forces)
    (; x, y, i, j, d2, d) = pair
    df = (1/d2)*(1/d)*(y - x)
    forces[i] += df
    forces[j] -= df
    return forces
end
map_pairwise!(update_forces!, system)
```

## Energy and forces

In this example, the potential energy and the forces are computed in a single
run, and a custom data structure is defined to store both values.

```julia
using CellListMap
using StaticArrays
# Define custom type
mutable struct EnergyAndForces
    energy::Float64
    forces::Vector{SVector{3,Float64}}
end
# Custom copy, reset and reducer functions
import CellListMap: copy_output, reset_output!, reducer
copy_output(x::EnergyAndForces) = EnergyAndForces(copy(x.energy), copy(x.forces))
function reset_output!(output::EnergyAndForces)
    output.energy = 0.0
    for i in eachindex(output.forces)
        output.forces[i] = SVector(0.0, 0.0, 0.0)
    end
    return output
end
function reducer(x::EnergyAndForces, y::EnergyAndForces)
    e_tot = x.energy + y.energy
    x.forces .+= y.forces
    return EnergyAndForces(e_tot, x.forces)
end
# Function that updates energy and forces for each pair
function energy_and_forces!(pair, output::EnergyAndForces)
    (; i, j, x, y, d2, d) = pair
    output.energy += 1/d
    df = (1/d2)*(1/d)*(y - x)
    output.forces[i] += df
    output.forces[j] -= df
    return output
end
# Initialize system
positions = rand(SVector{3,Float64},1000);
system = ParticleSystem(
    xpositions = positions,
    unitcell=[1.0,1.0,1.0],
    cutoff = 0.1,
    output = EnergyAndForces(0.0, similar(positions)),
    output_name = :energy_and_forces
)
# Compute energy and forces
map_pairwise(energy_and_forces!, system)
```

## Two sets of particles

In this example we illustrate the interface for the computation of properties
of two sets of particles, by computing the minimum distance between the two sets.

```julia
using CellListMap
using StaticArrays
# Custom structure to store the minimum distance pair
struct MinimumDistance
    i::Int
    j::Int
    d::Float64
end
# Function that updates the minimum distance found
function minimum_distance(pair, md)
    (; i, j, d) = pair
    if d < md.d
        md = MinimumDistance(i, j, d)
    end
    return md
end
# Define appropriate methods for copy, reset and reduce
import CellListMap: copy_output, reset_output!, reducer!
copy_output(md::MinimumDistance) = md
reset_output!(md::MinimumDistance) = MinimumDistance(0, 0, +Inf)
reducer!(md1::MinimumDistance, md2::MinimumDistance) = md1.d < md2.d ? md1 : md2
# Build system
xpositions = rand(SVector{3,Float64},100);
ypositions = rand(SVector{3,Float64},1000);
system = ParticleSystem(
       xpositions = xpositions,
       ypositions = ypositions,
       unitcell=[1.0,1.0,1.0],
       cutoff = 0.1,
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
# Compute the minimum distance
map_pairwise(minimum_distance, system)
```

In the above example, the function is used such that cell lists are constructed for both
sets. There are situations where this is not optimal, in particular:

1. When one of the sets if very small. In this case, constructing a cell list for the largest
   set becomes the bottleneck. Therefore, it is better to construct a cell list for the smallest
   set and loop over the particles of the largest set.
2. When one of the set is fixed and the second set is variable. In this case, it is better to
   construct the cell list for the fixed set only and loop over the variables of the variable set.

For dealing with these possibilities, an additional two-set interface is available, where one maps
the computation over an array of particles relative to a previously computed cell list. Complementing
the example above, we could compute the same minimum distance using:

```julia
# Construct the cell list system only for one of the sets: ypositions
ysystem = ParticleSystem(
       positions = ypositions,
       unitcell=[1.0,1.0,1.0],
       cutoff = 0.1,
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
# obtain the minimum distance between xpositions and the cell list in system
# Note the additional `xpositions` parameter in the call to map_pairwise.
map_pairwise(get_md, xpositions, ysystem)
```

Additionally, if the `xpositions` are updated, we can obtain compute the function relative to `ysystem` without
having to update the cell lists:

```julia-repl
julia> xpositions = rand(SVector{3,Float64},100);

julia> map_pairwise(minimum_distance, xpositions, ysystem)
MinimumDistance(67, 580, 0.008423693268450603)
```

while with the two-set cell list system one would need to update the cell lists for this new computation.

!!! compat
    The single-set cross-interaction was introduced in v0.10.0. It uses the method previously implemented
    for all cross-interactions.

## Particle simulation

In this example, a complete particle simulation is illustrated, with a simple potential.  This example can illustrate how particle positions and forces can be updated. Run this
simulation with:

```julia-repl
julia> system = init_system(N=200); # number of particles

julia> trajectory = simulate(system);

julia> animate(trajectory)
```

One important characteristic of this example is that the `system` is built outside the function that performs the simulation. This is done because the construction of the system is type-unstable (it is dimension, geometry and output-type dependent). Adding a function barrier avoids type-instabilities to propagate to the simulation causing possible performance problems.

```julia
using StaticArrays
using CellListMap
import CellListMap.wrap_relative_to
# Function that updates the forces, for potential of the form:
# if d < cutoff k*(d^2-cutoff^2)^2 else 0.0 with k = 10^6
function update_forces!(pair, forces, cutoff)
    (i, j, x, y, d2) = pair
    r = y - x
    dudr = 10^6 * 4 * r * (d2 - cutoff^2)
    forces[i] += dudr
    forces[j] -= dudr
    return forces
end
# Function that initializes the system: it is preferable to initialize
# the system outside the function that performs the simulation, because
# the system (data)type is defined on initialization. Initializing it outside
# the simulation function avoids possible type-instabilities.
function init_system(;N::Int=200)
    Vec2D = SVector{2,Float64}
    positions = rand(Vec2D, N)
    unitcell = [1.0, 1.0]
    cutoff = 0.1
    system = ParticleSystem(
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=similar(positions),
        output_name=:forces,
    )
    return system
end
function simulate(system=init_system(); nsteps::Int=100, isave=1)
    # initial velocities
    velocities = [ randn(eltype(system.positions)) for _ in 1:length(system.positions) ]
    dt = 1e-3
    trajectory = typeof(system.positions)[]
    for step in 1:nsteps
        # compute forces at this step
        map_pairwise!(
            (pair, forces) -> update_forces!(pair, forces, system.cutoff),
            system
        )
        # Update positions and velocities
        for i in eachindex(system.positions, system.forces)
            f = system.forces[i]
            x = system.positions[i]
            v = velocities[i]
            x = x + v * dt + (f / 2) * dt^2
            v = v + f * dt
            # wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)
            # !!! IMPORTANT: Update arrays of positions and velocities
            system.positions[i] = x
            velocities[i] = v
        end
        # Save step for printing
        if step % isave == 0
            push!(trajectory, copy(system.positions))
        end
    end
    return trajectory
end

using Plots
function animate(trajectory)
    anim = @animate for step in trajectory
        scatter(
            Tuple.(step),
            label=nothing,
            lims=(-0.5, 0.5),
            aspect_ratio=1,
            framestyle=:box,
        )
    end
    gif(anim, "simulation.gif", fps=10)
end
```
