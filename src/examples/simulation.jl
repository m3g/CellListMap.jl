using StaticArrays
using CellListMap
import CellListMap.wrap_relative_to
# Function that updates the forces, for potential of the form:
# if d < cutoff k*(d^2-cutoff^2)^2 else 0.0 with k = 10^6
function update_forces!(x, y, i, j, d2, forces, cutoff)
    r = y - x
    dudr = 10^6 * 4 * r * (d2 - cutoff^2)
    forces[i] += dudr
    forces[j] -= dudr
    return forces
end
# Function that initializes the system: it is preferrable to initialize
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
            (x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces,system.cutoff),
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