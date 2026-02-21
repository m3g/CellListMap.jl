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
pairwise!(minimum_distance, system)
```

