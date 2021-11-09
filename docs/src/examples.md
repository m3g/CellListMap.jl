# Examples

The full code of the examples described here is available at the [examples](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/) directory. 

## Mean difference of coordinates 

Computing the mean difference in `x` position between random particles. The closure is used to remove the indexes and the distance of the particles from the parameters of the input function, as they are not needed in this case.

```julia
using CellListMap

# System properties
N = 100_000
sides = [250,250,250]
cutoff = 10

# Particle positions
x = [ sides .* rand(3) for i in 1:N ]

# Initialize linked lists and box structures
box = Box(sides,cutoff)
cl = CellList(x,box)

# Function to be evaluated from positions 
f(x,y,sum_dx) = sum_dx + abs(x[1] - y[1])
normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

# Run calculation (0.0 is the initial value)
avg_dx = normalization * map_pairwise!(
    (x,y,i,j,d2,sum_dx) -> f(x,y,sum_dx), 0.0, box, cl 
)
```

The example above can be run with `CellListMap.Examples.average_displacement()` and is available in the
[average_displacement.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/average_displacement.jl) file.

## Histogram of distances

Computing the histogram of the distances between particles (considering the same particles as in the above example). Again,
we use a closure to remove the positions and indexes of the particles from the function arguments, because they are not needed. The distance, on the other side, is needed in this example:

```julia
# Function that accumulates the histogram of distances
function build_histogram!(d2,hist)
    d = sqrt(d2)
    ibin = floor(Int,d) + 1
    hist[ibin] += 1
    return hist
end;

# Initialize (and preallocate) the histogram
hist = zeros(Int,10);
normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

# Run calculation
hist = normalization * map_pairwise!(
    (x,y,i,j,d2,hist) -> build_histogram!(d2,hist),
    hist,box,cl
)

```

The example above can be run with `CellListMap.Examples.distance_histogram()` and is available in the
[distance_histogram.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/distance_histogram.jl) file.

## Gravitational potential

In this test we compute the "gravitational potential", assigning to each particle a different mass. In this case, the closure is used to pass the masses to the function that computes the potential.

```julia
# masses
const mass = rand(N)

# Function to be evaluated for each pair 
function potential(i,j,d2,mass,u)
    d = sqrt(d2)
    u = u - 9.8*mass[i]*mass[j]/d
    return u
end

# Run pairwise computation
u = map_pairwise!((x,y,i,j,d2,u) -> potential(i,j,d2,mass,u),0.0,box,cl)
```

The example above can be run with `CellListMap.Examples.gravitational_potential()` and is available in the
[gravitational_potential.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/gravitational_potential.jl) file.

## Gravitational force

In the following example, we update a force vector of for all particles.

```julia
# masses
const mass = rand(N)

# Function to be evaluated for each pair: update force vector
function calc_forces!(x,y,i,j,d2,mass,forces)
    G = 9.8*mass[i]*mass[j]/d2
    d = sqrt(d2)
    df = (G/d)*(x - y)
    forces[i] = forces[i] - df
    forces[j] = forces[j] + df
    return forces
end

# Initialize and preallocate forces
forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

# Run pairwise computation
forces = map_pairwise!(
    (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
    forces,box,cl
)

```

The example above can be run with `CellListMap.Examples.gravitational_force()` and is available in the
[gravitational_force.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/gravitational_force.jl) file.

## Nearest neighbour

Here we compute the indexes of the particles that satisfy the minimum distance between two sets of points, using the linked lists. The distance and the indexes are stored in a tuple, and a reducing method has to be defined for that tuple to run the calculation.  The function does not need the coordinates of the points, only their distance and indexes.

```julia
# Number of particles, sides and cutoff
N1=1_500
N2=1_500_000
sides = [250,250,250]
cutoff = 10.
box = Box(sides,cutoff)

# Particle positions
x = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:N1 ]
y = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:N2 ]

# Initialize auxiliary linked lists
cl = CellList(x,y,box)

# Function that keeps the minimum distance
f(i,j,d2,mind) = d2 < mind[3] ? (i,j,d2) : mind

# We have to define our own reduce function here
function reduce_mind(output,output_threaded)
    mind = output_threaded[1]
    for i in 2:Threads.nthreads()
        if output_threaded[i][3] < mind[3]
            mind = output_threaded[i]
        end
    end
    return mind
end

# Initial value
mind = ( 0, 0, +Inf )

# Run pairwise computation
mind = map_pairwise!( 
    (x,y,i,j,d2,mind) -> f(i,j,d2,mind),
    mind,box,cl;reduce=reduce_mind
)
```

The example above can be run with `CellListMap.Examples.nearest_neighbour()` and is available in the
[nearest_neighbour.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/nearest_neighbour.jl) file.

The example `CellListMap.Examples.nearest_neighbour_nopbc()` of [nearest\_neighbour\_nopbc.jl](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/nearest_neighbour_nopbc.jl) describes a similar problem but *without* periodic boundary conditions. Depending on the distribution of points and size it is a faster method than usual ball-tree methods. 

## Neighbour list

Obs: The package provides a `neighbourlist` function that implements this calculation, and can be used with:

```julia-repl
julia> x = [ rand(3) for i in 1:10_000 ];

julia> CellListMap.neighbourlist(x,0.05)
24778-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 62, 0.028481068525796384)
 ⋮
 (9954, 1749, 0.04887502372299809)
 (9974, 124, 0.040110356034451795)
```
or `CellListMap.neighbourlist(x,y,r)` for computing the lists of pairs of two sets closer than `r`.

The returning array contains tuples with the index of the particle in the first vector, the index of the particle in the second vector, and their distance.

The implementation of this function follows the principles below. 
 The empty `pairs` output array will be split in one vector for each thread, and reduced with a custom reduction function. 

```julia
# Function to be evaluated for each pair: push pair
function push_pair!(i,j,d2,pairs)
    d = sqrt(d2)
    push!(pairs,(i,j,d))
    return pairs
end

# Reduction function
function reduce_pairs(pairs,pairs_threaded)
    pairs = pairs_threaded[1]
    for i in 2:Threads.nthreads()
        append!(pairs,pairs_threaded[i])
    end
    return pairs
end

# Initialize output
pairs = Tuple{Int,Int,Float64}[]

# Run pairwise computation
pairs = map_pairwise!(
    (x,y,i,j,d2,pairs) -> push_pair!(i,j,d2,pairs),
    pairs,box,cl,
    reduce=reduce_pairs
)
```

The full example can be run with `CellListMap.Examples.neighbourlist()`, available in the file 
[neighbourlist.jl](https://github.com/m3g/CellListMap.jl/blob/src/examples/neighbourlist.jl).

## Units, automatic differentiation, etc.

The functions of CellListMap.jl support the propagation of generic (isbits) types, and thus units and thus automatic differentiation and the use of `Unitful`. A set of working examples can be found in the [generic_types.jl](https://github.com/m3g/CellListMap.jl/blob/src/examples/generic_types.jl) file.

### `Unitful` and units

We start illustrating the support for unit propagation. We need to define all involved quantities in the same units:

```julia-repl
julia> using Unitful, StaticArrays

julia> cutoff = 0.1u"nm" 
0.1 nm

julia> box = Box([1.0, 1.0, 1.0]u"nm",cutoff)
Box{OrthorhombicCell, 3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, Quantity{Float64, 𝐋^2, Unitful.FreeUnits{(nm^2,), 𝐋^2, nothing}}, 9}
  unit cell matrix = [ 1.0 nm, 0.0 nm, 0.0 nm; 0.0 nm, 1.0 nm, 0.0 nm; 0.0 nm, 0.0 nm, 1.0 nm ]
  cutoff = 0.1 nm
  number of computing cells on each dimension = [12, 12, 12]
  computing cell sizes = [0.1 nm, 0.1 nm, 0.1 nm] (lcell: 1)
  Total number of cells = 1728

julia> x = [ rand(typeof(cutoff),3) for _ in 1:1000 ];

julia> cl = CellList(x,box)
CellList{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}
  1000 real particles.
  626 cells with real particles.
  1694 particles in computing box, including images.
```

The corresponding mapping must take care of defining the result in the correct units associated to the expected output. For example, here we will compute just the sum of the squared distances between the particles within the cutoff. Thus, the expected output has the same units as the square of the dimensions: 

```julia-repl
julia> sum_sqr = zero(typeof(cutoff^2))
0.0 nm^2

julia> sum_sqr = map_pairwise!(
           (x,y,i,j,d2,sum_sqr) -> sum_sqr += d2,
           sum_sqr, box, cl
       )
12.983283925249138 nm^2
```

The performance penalty associated to propagating units is small. With units, we get:
```
julia> using BenchmarkTools

julia> @btime let x = $x
           cutoff = 0.1u"nm" 
           box = Box([1.0,1.0,1.0]u"nm",cutoff)
           cl = CellList(x,box)
           sum_sqr = zero(typeof(cutoff^2))
           map_pairwise!(
               (x,y,i,j,d2,sum_sqr) -> sum_sqr += d2,
               sum_sqr, box, cl
           )
       end
  1.828 ms (6547 allocations: 1.22 MiB)
12.983283925249138 nm^2
```

and the same problem without units runs in:

```julia-repl
julia> @btime let x = $([ SVector{3,Float64}(ustrip.(v)) for v in x ])
           cutoff = 0.1 
           box = Box([1.0,1.0,1.0],cutoff)
           cl = CellList(x,box)
           sum_sqr = 0.
           map_pairwise!(
               (x,y,i,j,d2,sum_sqr) -> sum_sqr += d2,
               sum_sqr, box, cl
           )
       end
  1.783 ms (6547 allocations: 1.22 MiB)
12.983283925249138
```

Auxiliary functions, like `CellListMap.neighbourlist`, propagate units correctly:

```julia-repl
julia> cutoff = 0.1u"nm"
0.1 nm

julia> box = Box([1.0, 1.0, 1.0]u"nm",cutoff);

julia> x = [ rand(typeof(cutoff),3) for _ in 1:1000 ];

julia> CellListMap.neighbourlist(x,cutoff)
1796-element Vector{Tuple{Int64, Int64, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}}:
 (1, 583, 0.06456224519583421 nm)
 (10, 216, 0.04958058924623024 nm)
 ⋮
 (934, 615, 0.08834318454969409 nm)
 (934, 692, 0.05002019032986014 nm)
```

### Automatic differentiation

Allowing automatic differentiation follows the same principles, meaning that we only need to allow the propagation of dual types through the computation by proper initialization of the input data.

The variables are each component of each vector, thus the easiest way to represent the points such that automatic differentiation packages understand is by creating a matrix:

```julia-repl
julia> x = rand(3,1000)
3×1000 Matrix{Float64}:
 0.186744  0.328719  0.874102  0.503535   …  0.328161  0.0895699  0.917338
 0.176157  0.972954  0.80729   0.624724      0.655268  0.470754   0.327578
 0.648482  0.537362  0.599624  0.0688776     0.92333   0.497984   0.208924
```

The key here is allow all the types of the parameters to follow the type propagation of the elements of `x` inside the differentiation routine. The function we define to compute the derivative is, then:

```julia-repl
julia> function sum_sqr(x,sides,cutoff)
           cutoff = eltype(x)(cutoff)
           sides = eltype(x).(sides)
           box = Box(sides,cutoff)
           cl = CellList(x,box)
           sum_sqr = zero(eltype(x))
           sum_sqr = map_pairwise!(
               (x,y,i,j,d2,sum_sqr) -> sum_sqr += d2,
               sum_sqr, box, cl
           )
           return sum_sqr
       end
sum_sqr (generic function with 1 method)
```

Note that we allow `cutoff`  and `sides`  to be converted to the same type of the input `x`  of the function. For a simple call to the function this is inconsequential:

```julia-repl
julia> cutoff = 0.1; sides = [1,1,1];

julia> sum_sqr(x,sides,cutoff)
12.897650398753228
```

but the conversion is required to allow the differentiation to take place:

```julia-repl
julia> ForwardDiff.gradient(x -> sum_sqr(x,sides,cutoff),x)
3×1000 Matrix{Float64}:
 -0.132567   0.029865  -0.101301  …   0.249267    0.0486424  -0.0400487
  0.122421   0.207495  -0.184366     -0.201648   -0.105031    0.218342
  0.0856502  0.288924   0.122445     -0.0147022  -0.103314   -0.0862264
```

### Measurements

Propagating uncertainties through the `Measurements`  and other similar packages requires a different strategy, because within `CellListMap` only `isbits` types can be used, which is not the case of the type `Measurement` type. 

In cases like this, it is better to bypass all the internals of `CellListMap`  and provide the data to the function that computes pairwise properties directly as a closure. For example:

A vector of particles with uncertainties in their coordinates can be created with: 
```julia-repl
julia> using StaticArrays 

julia> x_input = [ SVector{3}(measurement(rand(),0.01*rand()) for i in 1:3) for j in 1:1000 ]
1000-element Vector{SVector{3, Measurement{Float64}}}:
 [0.1658 ± 0.003, 0.9951 ± 0.0054, 0.5067 ± 0.0035]
 [0.2295 ± 0.0074, 0.2987 ± 0.0021, 0.42828 ± 0.00099]
 ⋮
 [0.1362 ± 0.0034, 0.2219 ± 0.0048, 0.2119 ± 0.0072]
 [0.2521 ± 0.0038, 0.4988 ± 0.00013, 0.856046 ± 4.3e-5]
```

The variables within the `CellListMap` functions will be stripped from the uncertainties. We do:

```julia-repl
julia> cutoff = 0.1; box = Box([1,1,1],cutoff);

julia> x_strip = [ getproperty.(v,:val) for v in x_input ]
1000-element Vector{SVector{3, Float64}}:
 [0.08441931492362276, 0.9911530546181084, 0.07408559584648788]
 [0.12084764467339837, 0.8284551316333133, 0.9021906852432111]
 ⋮
 [0.2418752113326077, 0.4429225751775432, 0.13576355747772784]
 [0.24440380524702654, 0.07148275176890073, 0.26722687487212315]
 ```

The cell list is built with the stripped values:

```julia-repl
julia> cl = CellList(x_strip,box)
CellList{3, Float64}
  1000 real particles.
  637 cells with real particles.
  1695 particles in computing box, including images.
```

The result is initialized with the proper type,

```
julia> sum_sqr = measurement(0.,0.)
0.0 ± 0.0
```

and the mapping is performed with the stripped coordinates, but passing the values with uncertainties to the mapped function, which will perform the computation on the pairs with those values:

```
julia> using LinearAlgebra: norm_sqr

julia> sum_sqr = map_pairwise!(
           (xᵢ,xⱼ,i,j,d2,sum_sqr) -> begin
               x1 = x_input[i]
               x2 = CellListMap.wrap_relative_to(x_input[j],x1,box)
               sum_sqr += norm_sqr(x2-x1)
               return sum_sqr
           end, 
           sum_sqr, box, cl
       )
13.14 ± 0.061
```

In the function above, the `xᵢ` and `xⱼ` coordinates, which correspond to the coordinates in `x_input[i]` and `x_input[j]`, but already wrapped relative to each other, are ignored, because they don't carry the uncertainties. We use only the indexes `i` and `j` to recompute the relative position of the particles according to the periodic boundary conditions (using the `CellListMap.wrap_relative_to` function) and their (squared) distance. Since the `x_input`  array carries the uncertainties, the computation of `sum_sqr` will propagate them.   

!!! note
    All these computations should be performed inside the scope of a function for optimal performance. The examples here can be followed by copying and pasting the code into the REPL, but this is not the recommended practice for critical code. The strategy of bypassing the internal computations of `CellListMap` may be useful for improving performance even if the previous and simpler method is possible. 































