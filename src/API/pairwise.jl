"""
    pairwise!(
        f::Function, system::AbstractParticleSystem; 
        show_progress=true, reset=true,
    )

Function that maps the `f` function into all pairs of particles of
`system` that are found to be within the `cutoff`. 

The function `f` receives a `NeighborPair` struct and the output:
```
function f(pair, output)
    # pair.i, pair.j: indices of the particles
    # pair.x, pair.y: coordinates (minimum-image adjusted)
    # pair.d: distance between particles
    # pair.d2: squared distance
    # update output
    return output
end
```

Thread-safety is taken care automatically in parallel executions.

`pairwise` is an alias to `pairwise!` for syntax consistency
when the `output` variable is immutable.

If `reset` is set to `false`, the value of `system.output` will not be
set to `zero(typeof(system.output))` before the new accumulation.

# Example

In this example we compute the sum of `1/(1+d)` where `d` is the
distance between particles of a set, for `d < cutoff`. 

```julia-repl
julia> sys = ParticleSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0
           );

julia> pairwise!((pair, output) -> output += 1 / (1 + pair.d), sys)
1870.0274887950268
```

"""
function pairwise!(
        f::F,
        sys::AbstractParticleSystem;
        show_progress::Bool = false,
        reset::Bool = true,
    ) where {F <: Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded; reset)
    UpdateParticleSystem!(sys)
    sys.output = _pairwise!(
        f, sys.output, sys._box, sys._cell_list;
        output_threaded = sys._output_threaded,
        parallel = sys.parallel,
        show_progress = show_progress
    )
    return sys.output
end

#
# Cross-computations when only one cell list was computed
#
"""
    pairwise!(f::Function, x::AbstractVector{<:AbstractVector}, sys::ParticleSystem1; kargs...)
    pairwise!(f::Function, x::AbstractMatrix, sys::ParticleSystem1; kargs...)

Evaluate function f for pairs in two independent sets of particles, where the `sys::ParticleSystem1` object
contains the previously computed cell lists of one set of particles, and the second set is given by the
array of positions `x`.

This function can be advantageous over computing the interactions with `CellListPair`, because here the
cell lists are only computed for one set. This may be advantageous in two situations:

1. The second set of particles is not changing, and the first set is changing. Thus, the cell lists
   of the second set can be computed only once.
2. One of the sets is much smaller than the other. In this case, computing the cell lists of the largest
   set might be too expensive. Construct the `ParticleSystem` object for the smallest set, and use this
   function to compute the interactions with the largest set.

## Keyword arguments:

- `show_progress::Bool=false`: Show progress bar.
- `reset::Bool=true`: If set to `false` the value of `sys.output` will not be set to `zero(typeof(sys.output)`,
   and the result will be accumulated

## Example

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 1000);

julia> sys = ParticleSystem(positions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0);

julia> y = rand(SVector{3,Float64}, 100);

julia> pairwise!((pair, output) -> output + pair.d, y, sys) # Compute the sum of the distances of x and y
31.121496300032163

julia> z = rand(SVector{3,Float64}, 200);

julia> pairwise!((pair, output) -> output + pair.d, z, sys) # Compute the sum of the distances x and z
63.57860511891242
```

"""
function pairwise!(
        f::F, x::AbstractVecOrMat, sys::ParticleSystem1;
        show_progress::Bool = false, reset::Bool = true,
) where {F <: Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded; reset)
    UpdateParticleSystem!(sys)
    sys.output = _pairwise!(
        f, sys.output, sys._box, x, sys._cell_list;
        output_threaded = sys._output_threaded,
        parallel = sys.parallel, show_progress,
    )
    return sys.output
end
