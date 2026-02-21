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

