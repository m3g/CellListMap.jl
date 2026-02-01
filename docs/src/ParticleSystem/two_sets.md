# Two sets of particles

If the computation involves two sets of particles, a similar interface is available.
The only difference is that the coordinates of the two sets must be provided to
the `ParticleSystem` constructor as the `xpositions` and `ypositions` arrays.

We will illustrate this interface by computing the minimum distance between two
sets of particles, which allows us to showcase further the definition of custom
type interfaces:

First, we define a variable type that will carry the indexes and
the distance of the closest pair of particles:
```julia-repl
julia> struct MinimumDistance
           i::Int
           j::Int
           d::Float64
       end
```

The function that, given two particles, retains the minimum distance, is:
```julia-repl
julia> function minimum_distance(i, j, d2, md)
           d = sqrt(d2)
           if d < md.d
               md = MinimumDistance(i, j, d)
           end
           return md
       end
minimum_distance (generic function with 1 method)
```

We overload copy, reset, and reduce functions, accordingly:
```julia-repl
julia> import CellListMap: copy_output, reset_output!, reducer!

julia> copy_output(md::MinimumDistance) = md
copy_output (generic function with 5 methods)

julia> reset_output!(md::MinimumDistance) = MinimumDistance(0, 0, +Inf)
reset_output! (generic function with 5 methods)

julia> reducer!(md1::MinimumDistance, md2::MinimumDistance) = md1.d < md2.d ? md1 : md2
reducer! (generic function with 2 methods)
```
Note that since `MinimumDistance` is immutable, copying it is the same as returning the value.
Also, resetting the minimum distance consists of setting its `d` field to `+Inf`. And, finally,
reducing the threaded distances consists of keeping the pair with the shortest distance.

Next, we build the system

```julia-repl
julia> xpositions = rand(SVector{3,Float64},1000);

julia> ypositions = rand(SVector{3,Float64},1000);

julia> system = ParticleSystem(
           xpositions = xpositions,
           ypositions = ypositions,
           unitcell=[1.0,1.0,1.0],
           cutoff = 0.1,
           output = MinimumDistance(0,0,+Inf),
           output_name = :minimum_distance,
        )
```

And finally we can obtain the minimum distance between the sets:

```julia-repl
julia> map_pairwise((x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), system)
MinimumDistance(276, 617, 0.006009804808785543)
```

## Cross-interactions with a single cell list

The two-set interface above builds cell lists for both sets. There are situations where this is not optimal:

1. **One set is fixed, the other changes.** If `ypositions` never changes but `xpositions` is updated each step,
   rebuilding a cell list for `ypositions` every time is wasteful. Build it once and reuse.
2. **One set is much smaller.** Building a cell list for a very large set can be expensive. Build the cell list
   only for the smaller set and iterate over the larger set directly.

For these cases, construct a `ParticleSystem` for the *reference* set (the one whose cell list you want to keep),
then pass the *other* set as a plain array to `map_pairwise`:

```julia-repl
julia> # Build the system only for ypositions (the reference set)

julia> ysystem = ParticleSystem(
           positions = ypositions,
           unitcell = [1.0, 1.0, 1.0],
           cutoff = 0.1,
           output = MinimumDistance(0, 0, +Inf),
           output_name = :minimum_distance,
       )

julia> # Compute interactions between xpositions and the cell list in ysystem.
       # Note: xpositions is passed as a plain array, before the system argument.

julia> map_pairwise((x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), xpositions, ysystem)
MinimumDistance(67, 580, 0.008423693268450603)
```

When `xpositions` changes, the cell list in `ysystem` does not need to be rebuilt:

```julia-repl
julia> xpositions = rand(SVector{3,Float64}, 100);  # new positions

julia> map_pairwise((x,y,i,j,d2,md) -> minimum_distance(i,j,d2,md), xpositions, ysystem; update_lists=false)
MinimumDistance(42, 310, 0.005182736498125341)
```

The `update_lists=false` keyword skips updating the cell list of `ysystem`, since only `xpositions` changed.
If `ysystem.positions` itself changed, use `update_lists=true` (the default).

!!! compat
    The single-set cross-interaction interface (`map_pairwise(f, x, system)`) was introduced in v0.10.0.

## Benchmarking of cross-interaction alternatives

With the following functions we will benchmark the performance of the two alternatives for computing
cross-set interactions, **including** the time required to build the cell lists (the initialization
of the `ParticleSystem`) objects:

```julia
# First alternative: compute cell lists for the two sets
function two_set_celllist(xpositions, ypositions)
   system = ParticleSystem(
       xpositions = xpositions,
       ypositions = ypositions,
       unitcell=[1.0,1.0,1.0],
       cutoff = 0.1,
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
    return map_pairwise(get_md, system)
end
# Second alternative: compute cell lists for one set
function one_set_celllist(xpositions, ypositions)
   system = ParticleSystem(
       positions = ypositions,
       unitcell=[1.0,1.0,1.0],
       cutoff = 0.1,
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
)
    return map_pairwise(get_md, xpositions, system)
end
```

If one of the sets is small, the one-set alternative is clearly faster, if we
construct the cell lists for the smaller set:

```julia-repl
julia> using BenchmarkTools

julia> xpositions = rand(SVector{3,Float64}, 10^6);

julia> ypositions = rand(SVector{3,Float64}, 100);

julia> @btime one_set_celllist($xpositions, $ypositions) samples=1 evals=1
  25.165 ms (1531 allocations: 575.72 KiB)
MinimumDistance(65937, 63, 0.00044803040276614203)

julia> @btime two_set_celllist($xpositions, $ypositions) samples=1 evals=1
  207.129 ms (154794 allocations: 478.00 MiB)
MinimumDistance(65937, 63, 0.00044803040276614203)
```

For much larger system, though, the computation of the cell lists become less relevant and the first alternative
might be the most favorable, even including the cell lists updates:

```julia-repl
julia> @btime one_set_celllist($xpositions, $ypositions) samples=1 evals=1
  12.196 s (153327 allocations: 478.02 MiB)
MinimumDistance(627930, 889247, 7.59096139675071e-5)

julia> @btime two_set_celllist($xpositions, $ypositions) samples=1 evals=1
  2.887 s (306416 allocations: 952.00 MiB)
MinimumDistance(627930, 889247, 7.59096139675071e-5)
```

This performance advantage of the two-set cell lists arises because more interactions can be skipped
by [precomputing properties of the cells involved](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20563).
On the other side, when the lists are available
for only one set, the loop over all the particles of the second set is mandatory. Since this loop
is fast, it is favorable over the construction of the cell lists for smaller sets.
