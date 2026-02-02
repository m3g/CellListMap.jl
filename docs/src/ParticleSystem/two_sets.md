# Two sets of particles

If the computation involves two sets of particles, a similar interface is available.
The only difference is that the coordinates of the two sets must be provided to
the `ParticleSystem` constructor as the `xpositions` and `ypositions` arrays.

## Minimum-distance example

We will illustrate this interface by computing the minimum distance between two
sets of particles, which allows us to showcase further the definition of custom
type interfaces:

```@example mind
using CellListMap, StaticArrays
```

First, we define a variable type that will carry the indexes and
the distance of the closest pair of particles:

```@example mind
struct MinimumDistance
    i::Int
    j::Int
    d::Float64
end
```

The function that, given two particles, retains the minimum distance, is:
```@example mind
function minimum_distance(pair, md)
    (; i, j, d) = pair
    if d < md.d
        md = MinimumDistance(i, j, d)
    end
    return md
end
```

We overload `copy`, `reset`, and `reducer` functions, accordingly:
```@example mind
import CellListMap: copy_output, reset_output!, reducer!
copy_output(md::MinimumDistance) = md
```
```@example mind
reset_output!(md::MinimumDistance) = MinimumDistance(0, 0, +Inf)
```
```@example mind
reducer!(md1::MinimumDistance, md2::MinimumDistance) = md1.d < md2.d ? md1 : md2
```

Note that since `MinimumDistance` is immutable, copying it is the same as returning the value.
Also, resetting the minimum distance consists of setting its `d` field to `+Inf`. And, finally,
reducing the threaded distances consists of keeping the pair with the shortest distance.

## Cross-interactions with two cell lists

Next, we build the system. Here we choose to provide both sets of particles to the `ParticleSystem` constructor, which means that cell lists will be built for both sets:

```@example mind
xpositions = rand(SVector{3,Float64},1000);
ypositions = rand(SVector{3,Float64},1000);
system = ParticleSystem(
   xpositions = xpositions,
   ypositions = ypositions,
   unitcell=[1.0,1.0,1.0],
   cutoff = 0.1,
   output = MinimumDistance(0,0,+Inf),
   output_name = :minimum_distance,
)
```

And finally we can obtain the minimum distance between the sets:

```@example mind
map_pairwise!(minimum_distance, system)
```

## Cross-interactions with a single cell list

The two-set interface above builds cell lists for both sets. There are situations where this is not optimal:

1. **One set is fixed, the other changes.** If `ypositions` never changes but `xpositions` is updated each step, rebuilding a cell list for `ypositions` every time is wasteful. Build it once and reuse.
2. **One set is much smaller.** Building a cell list for a very large set can be expensive. Build the cell list only for the smaller set and iterate over the larger set directly.

For these cases, construct a `ParticleSystem` for the *reference* set (the one whose cell list you want to keep), then pass the *other* set as a plain array to `map_pairwise`:

Build the system only for ypositions (the reference set):

```@example mind
ysystem = ParticleSystem(
    positions = ypositions,
    unitcell = [1.0, 1.0, 1.0],
    cutoff = 0.1,
    output = MinimumDistance(0, 0, +Inf),
    output_name = :minimum_distance,
)
```

Compute interactions between xpositions and the cell list in ysystem.
Note: `xpositions` is passed as a plain array, before the system argument.

```@example mind
map_pairwise!(minimum_distance, xpositions, ysystem)
```

When `xpositions` changes, the cell list in `ysystem` does not need to be rebuilt:

```@example mind
xpositions = rand(SVector{3,Float64}, 100);  # new positions
map_pairwise!(minimum_distance, xpositions, ysystem; update_lists=false)
```

The `update_lists=false` keyword skips updating the cell list of `ysystem`, since only `xpositions` changed. If `ysystem.positions` itself changed, use `update_lists=true` (the default).

!!! compat
    The single-set cross-interaction interface (`map_pairwise!(f, x, system)`) was introduced in v0.10.0.

## Benchmarking of cross-interaction alternatives

With the following functions we will benchmark the performance of the two alternatives for computing
cross-set interactions, **including** the time required to build the cell lists (the initialization
of the `ParticleSystem`) objects:

### First alternative: compute cell lists for the two sets

```@example mind
function two_set_celllist(xpositions, ypositions)
   system = ParticleSystem(
       xpositions = xpositions,
       ypositions = ypositions,
       unitcell=[1.0,1.0,1.0],
       cutoff = 0.1,
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
    )
    return map_pairwise!(minimum_distance, system)
end
```

### Second alternative: compute cell lists for one set

```@example mind
function one_set_celllist(xpositions, ypositions)
   system = ParticleSystem(
       positions = ypositions,
       unitcell=[1.0,1.0,1.0],
       cutoff = 0.1,
       output = MinimumDistance(0,0,+Inf),
       output_name = :minimum_distance,
    )
    return map_pairwise!(minimum_distance, xpositions, system)
end
```

If one of the sets is small, the one-set alternative is clearly faster, if we
construct the cell lists for the smaller set:

```@example mind
using BenchmarkTools
xpositions = rand(SVector{3,Float64}, 10^5);
ypositions = rand(SVector{3,Float64}, 100);
@benchmark one_set_celllist($xpositions, $ypositions) samples=1 evals=1
```

```@example mind
@benchmark two_set_celllist($xpositions, $ypositions) samples=1 evals=1
```

If both sets are large, though, the computation of the cell lists become less relevant and the first alternative might be the most favorable, even including the cell lists updates:

```@example mind
xpositions = rand(SVector{3,Float64}, 10^5);
ypositions = rand(SVector{3,Float64}, 10^5);
@benchmark one_set_celllist($xpositions, $ypositions) samples=1 evals=1
```

```@example mind
@benchmark two_set_celllist($xpositions, $ypositions) samples=1 evals=1
```

This performance advantage of the two-set cell lists arises because more interactions can be skipped by [precomputing properties of the cells involved](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20563). On the other side, when the lists are available
for only one set, the loop over all the particles of the second set is mandatory. Since this loop
is fast, it is favorable over the construction of the cell lists for smaller sets.
