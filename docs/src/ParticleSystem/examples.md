# Complete examples

This section contains complete, self-contained example codes that can be copied
and run directly.

## Simple energy computation

In this example, a simple potential energy defined as the sum of the
inverse of the distance between the particles is computed.

```@example
using CellListMap
using StaticArrays
system = ParticleSystem(
    xpositions = rand(SVector{3,Float64},1000),
    unitcell=[1.0,1.0,1.0],
    cutoff = 0.1,
    output = 0.0,
    output_name = :energy
)
pairwise!((pair, energy) -> energy += 1 / pair.d, system)
println(" Energy: $(system.energy)")
```

## Force computation

Here we compute the force vector associated to the potential energy
function of the previous example.

```@example
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
pairwise!(update_forces!, system)
println(""" 
    Force on particle 1: $(system.forces[1]))
    Force on particle 2: $(system.forces[1]))
""")
```

## Energy and forces

In this example, the potential energy and the forces are computed in a single
run, and a custom data structure is defined to store both values.

```@example 
using CellListMap
using StaticArrays
# Define custom type to store energy and force vector
mutable struct EnergyAndForces
    energy::Float64
    forces::Vector{SVector{3,Float64}}
end
# Custom copy, reset and reducer functions
import CellListMap: copy_output, reset_output!, reducer
copy_output(x::EnergyAndForces) = EnergyAndForces(copy(x.energy), copy(x.forces))
function reset_output!(x::EnergyAndForces)
    x.energy = 0.0
    fill!(x.forces, SVector(0.0, 0.0, 0.0))
    return x
end
function reducer(x::EnergyAndForces, y::EnergyAndForces)
    x.energy += y.energy
    x.forces .+= y.forces
    return x
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
pairwise!(energy_and_forces!, system)
# Print some results
println("""
    Energy: $(system.energy_and_forces.energy)
    Force on particle 1: $(system.energy_and_forces.forces[1])
""")
```

## Two sets of particles

In this example we illustrate the interface for the computation of properties
of two sets of particles, by computing the minimum distance between the two sets.

```@example
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
    md = d < md.d ? MinimumDistance(i, j, d) : md
    return md
end
# Define appropriate methods for copy, reset and reduce
import CellListMap: copy_output, reset_output!, reducer!
copy_output(md::MinimumDistance) = md
reset_output!(::MinimumDistance) = MinimumDistance(0, 0, +Inf)
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
pairwise!(minimum_distance, system)
```

## Particle simulation

In this example, a complete particle simulation is illustrated, with a simple potential.  This example can illustrate how particle positions and forces can be updated. Run this
simulation with:

```julia-repl
julia> trajectory = simulate(200) # 200 particles

julia> animate(trajectory)
```

Forces are computed similarly to what is done in the [Force computation](@ref) example. Additionally, we use
the public (but not exported) `CellListMap.wrap_relative_to` function, to keep the saved coordinates within
the minimal periodic coordinates relative to the origin, to produce a nice animation at the end.

```@example ex_simulation
using StaticArrays
using CellListMap
import CellListMap.wrap_relative_to
# Function that updates the forces, for potential of the form:
# if d < cutoff k*(d^2-cutoff^2)^2 else 0.0 with k = 10^6
function update_forces!(pair, forces, cutoff)
    (;i, j, x, y, d2) = pair
    r = y - x
    dudr = 10^6 * 4 * r * (d2 - cutoff^2)
    forces[i] += dudr
    forces[j] -= dudr
    return forces
end
function simulate(N; nsteps::Int=100, isave=1)
    # initial velocities
    velocities = randn(SVector{2,Float64}, N)
    # Initial positions
    positions = rand(SVector{2,Float64}, N)
    # force array
    forces = zeros(SVector{2,Float64}, N)
    # initialize ParticleSystem
    system = ParticleSystem(
        positions=positions,
        cutoff=0.1,
        unitcell=[1.0, 1.0],
        output=forces,
        output_name=:forces,
    )
    dt = 1e-3
    trajectory = typeof(positions)[]
    for step in 1:nsteps
        # compute forces at this step
        pairwise!(
            (pair, forces) -> update_forces!(pair, forces, system.cutoff),
            system
        )
        # Update positions and velocities
        for i in eachindex(positions, velocities, system.forces)
            f = system.forces[i]
            x = positions[i]
            v = velocities[i]
            x = x + v * dt + (f / 2) * dt^2
            v = v + f * dt
            # wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)
            # update positions and velocities
            positions[i] = x
            velocities[i] = v
        end
        # Update ParticleSystem positions
        update!(system; positions=positions)
        # Save step for printing
        if step % isave == 0
            push!(trajectory, copy(positions))
        end
    end
    return trajectory
end
```

The following function will create an animation from the resulting trajectory of the simulation:

```@example ex_simulation
using Plots
function animate(trajectory)
    anim = @animate for step in trajectory
        scatter(
            Tuple.(step),
            label=nothing,
            lims=(-0.5, 0.5),
            aspect_ratio=1,
            framestyle=:box,
            size=(400,400),
            ticks=nothing,
        )
    end
    gif(anim, "simulation.gif", fps=10)
end
```

Now, running the simulation and creating an animation:

```@example ex_simulation
trajectory = simulate(200) # 200 particles
animate(trajectory)
```
