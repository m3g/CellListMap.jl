module CellListMap

using DocStringExtensions
using ProgressMeter
using Parameters
using StaticArrays
using Setfield
using Base.Threads
using LinearAlgebra: dot, norm, norm_sqr, cross

export Box
export CellList, UpdateCellList!
export map_pairwise!
export limits
export TriclinicCell
export OrthorhombicCell

include("./StructTypes.jl")
include("./CellOperations.jl")
include("./CoreComputing.jl")

"""

```
map_pairwise!(
  f::Function,
  output,
  box::Box,
  cl::CellList
  ;parallel::Bool=true,
  show_progress::Bool=false
)
```

This function will run over every pair of particles which are closer than 
`box.cutoff` and compute the Euclidean distance between the particles, 
considering the periodic boundary conditions given in the `Box` structure. 
If the distance is smaller than the cutoff, a function `f` of the 
coordinates of the two particles will be computed. 

The function `f` receives six arguments as input: 
```
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the 
second particle, the index of the first particle, the index of the second 
particle, the squared distance between them, and the `output` variable. 
It has also to return the same `output` variable. Thus, `f` may or not 
mutate `output`, but in either case it must return it. With that, it is 
possible to compute an average property of the distance of the particles 
or, for example, build a histogram. The squared distance `d2` is computed 
internally for comparison with the 
`cutoff`, and is passed to the `f` because many times it is used for the 
desired computation. 

## Example

Computing the mean absolute difference in `x` position between random particles, 
remembering the number of pairs of `n` particles is `n(n-1)/2`. The function does 
not use the indices or the distance, such that we remove them from the parameters 
by using a closure.

```julia-repl
julia> n = 100_000;

julia> box = Box([250,250,250],10);

julia> x = [ SVector{3,Float64}(sides .* rand(3)) for i in 1:n ];

julia> cl = CellList(x,box);

julia> f(x,y,sum_dx) = sum_dx + abs(x[1] - y[1])

julia> normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

julia> avg_dx = normalization * map_parwise!((x,y,i,j,d2,sum_dx) -> f(x,y,sum_dx), 0.0, box, cl)

```

"""
function map_pairwise!(f::F, output, box::Box, cl; 
    # Parallelization options
    parallel::Bool=true,
    output_threaded=(parallel ? [ deepcopy(output) for i in 1:nthreads() ] : nothing),
    reduce::Function=reduce,
    show_progress::Bool=false,
) where {F} # Needed for specialization for this function (avoids some allocations)
    if parallel && nthreads() > 1
        output = map_pairwise_parallel!(f,output,box,cl;
            output_threaded=output_threaded,
            reduce=reduce,
            show_progress=show_progress
        )
    else
        output = map_pairwise_serial!(f,output,box,cl,show_progress=show_progress)
    end
    return output
end

"""

```
map_pairwise!(f::Function,output,box::Box,cl::CellListPair)
```

The same but to evaluate some function between pairs of the particles of the vectors.

"""
function map_pairwise!(f::F1, output, box::Box, cl::CellListPair;
    # Parallelization options
    parallel::Bool=true,
    output_threaded=(parallel ? [ deepcopy(output) for i in 1:nthreads() ] : nothing),
    reduce::F2=reduce,
    show_progress::Bool=false
) where {F1,F2} # Needed for specialization for this function (avoids some allocations) 
    if cl.swap 
        fswap(x,y,i,j,d2,output) = f(y,x,j,i,d2,output) 
    else
        fswap = f
    end
    if parallel && nthreads() > 1
        output = map_pairwise_parallel!(fswap,output,box,cl;
            output_threaded=output_threaded,
            reduce=reduce,
            show_progress=show_progress
        )
    else
        output = map_pairwise_serial!(fswap,output,box,cl,show_progress=show_progress)
    end
    return output
end

# Utils
include("./utils.jl")

#
# Test and example functions
#
include("./examples.jl")
include("./testing.jl")

end # module



