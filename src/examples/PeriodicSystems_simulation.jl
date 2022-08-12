using StaticArrays
using CellListMap.PeriodicSystems
import CellListMap.wrap_relative_to
# Function that updates the forces (potential of the form (d^2-radius)^2
function update_forces!(x, y, i, j, d2, forces, radius)
    dx = x - y
    if d2 < radius^2
        dfdx = 4 * dx * (d2 - radius)
        forces[i] -= dfdx
        forces[j] += dfdx
    end
    return forces
end
function simulate(; N::Int=1000, nsteps::Int=100, isave=1)
    Vec2D = SVector{2,Float64}
    positions = rand(Vec2D, N)
    unitcell = [1.0, 1.0]
    system = PeriodicSystem(
        positions=positions,
        cutoff=0.1,
        unitcell=unitcell,
        output=similar(positions),
        output_name=:forces,
    )
    velocities = [-1.0 .+ 2.0 * randn(Vec2D) for _ in 1:N]
    dt = 1e-3
    radius = 0.5
    trajectory = typeof(positions)[]
    for step in 1:nsteps
        # compute forces at this step
        map_pairwise!(
            (x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces,radius),
            system
        )
        # Update positions and velocities
        for i in eachindex(system.forces)
            f = system.forces[i]
            x = system.positions[i]
            v = velocities[i]
            x = x + v * dt + (f / 2) * dt^2
            v = v + f * dt
            # wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(0.0, 0.0), unitcell)
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