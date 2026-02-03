#=
    reduce(output, output_threaded)

Most common reduction function, which sums the elements of the output. 
Here, `output_threaded` is a vector containing `nbatches(cl)` copies of
the `output` variable (a scalar or an array). Custom reduction functions 
must replace this one if the reduction operation is not a simple sum. 
The `output_threaded` array is, by default, created automatically by copying
the given `output` variable `nbatches(cl)` times. 

## Examples

Scalar reduction: 

```julia-repl
julia> output = 0.; output_threaded = [ 1, 2 ];

julia> CellListMap.reduce(output,output_threaded)
3
```

Array reduction:

```julia-repl
julia> output = [0,0]; output_threaded = [ [1,1], [2,2] ];

julia> CellListMap.reduce(output,output_threaded)
2-element Vector{Int64}:
 3
 3

julia> output
2-element Vector{Int64}:
 3
 3
```

=#
reduce(output::T, output_threaded::Vector{T}) where {T} = sum(output_threaded)
function reduce(output::T, output_threaded::Vector{T}) where {T <: AbstractArray}
    for ibatch in eachindex(output_threaded)
        @. output += output_threaded[ibatch]
    end
    return output
end
function reduce(output, output_threaded)
    T = typeof(output)
    throw(
        ArgumentError(
            """\n
            No method matching reduce(::$(typeof(output)),::$(typeof(output_threaded)))

            Please provide a method that appropriately reduces a `Vector{$T}`, with
            the signature:

            ```
            custom_reduce(output::$T, output_threaded::Vector{$T})
            ```

            The reduction function **must** return the `output` variable, even 
            if it is mutable.  

            See: https://m3g.github.io/CellListMap.jl/stable/parallelization/#Custom-reduction-functions

            """
        )
    )
end

#
# Auxiliary functions to control the exhibition of the progress meter
#
_next!(::Nothing) = nothing
_next!(p) = next!(p)

#
# Functions necessary for the projection/partition scheme
#
#=
    partition!(by, x::AbstractVector)

# Extended help

Function that reorders `x` vector by putting in the first positions the
elements with values satisfying `by(el)`. Returns the number of elements
that satisfy the condition.

=#
function partition!(by, x::AbstractVector)
    iswap = 1
    @inbounds for i in eachindex(x)
        if by(x[i])
            if iswap != i
                x[iswap], x[i] = x[i], x[iswap]
            end
            iswap += 1
        end
    end
    return iswap - 1
end

#=
    project_particles!(projected_particles,cellⱼ,cellᵢ,Δc,Δc_norm,box)

# Extended help

Projects all particles of the cell `cellⱼ` into unnitary vector `Δc` with direction 
connecting the centers of `cellⱼ` and `cellᵢ`. Modifies `projected_particles`, and 
returns a view of `projected particles, where only the particles for which
the projection on the direction of the cell centers still allows the particle
to be within the cutoff distance of any point of the other cell.

=#
function project_particles!(
        projected_particles, cellⱼ, cellᵢ,
        Δc, Δc_norm, box::Box{UnitCellType, N}
    ) where {UnitCellType, N}
    if box.lcell == 1
        margin = box.cutoff + Δc_norm / 2 # half of the distance between centers
    else
        margin = box.cutoff * (1 + sqrt(N) / 2) # half of the diagonal of the cutoff-box
    end
    iproj = 0
    @inbounds for j in 1:cellⱼ.n_particles
        pⱼ = cellⱼ.particles[j]
        xproj = dot(pⱼ.coordinates - cellᵢ.center, Δc)
        if abs(xproj) <= margin
            iproj += 1
            projected_particles[iproj] = ProjectedParticle(pⱼ.index, xproj, pⱼ.coordinates, pⱼ.real)
        end
    end
    pp = @view(projected_particles[1:iproj])
    return pp
end
