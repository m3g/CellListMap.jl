module CellListMap

using DocStringExtensions: TYPEDEF, TYPEDFIELDS 
using TestItems: @testitem
using ProgressMeter: Progress, next!
using Parameters: @unpack, @with_kw
using StaticArrays: SVector, SMatrix, @SVector, @SMatrix, MVector, MMatrix
using Setfield: @set!
using LinearAlgebra: cross, diagm, I
using Base.Threads: nthreads, @spawn 

export Box
export CellList, UpdateCellList!
export map_pairwise!, map_pairwise
export limits
export TriclinicCell
export OrthorhombicCell
export NonPeriodicCell
export unitcelltype
export nbatches

# Flag for internal function doc entries
const INTERNAL = "Internal function or structure - interface may change."

# Testing file
const argon_pdb_file = joinpath("$(@__DIR__ )/../test/gromacs/argon/cubic.pdb")


include("./linearalgebra.jl")
include("./show.jl")
include("./Box.jl")
include("./CellLists.jl")
include("./PeriodicSystems.jl")
include("./CellOperations.jl")
include("./CoreComputing.jl")

"""
    map_pairwise!(
        f::Function,
        output,
        box::Box,
        cl::CellList
        ;parallel::Bool=true,
        show_progress::Bool=false
    )

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
function map_pairwise!(f::F, output, box::Box, cl::CellList; 
    # Parallelization options
    parallel::Bool=true,
    output_threaded=nothing,
    reduce::Function=reduce,
    show_progress::Bool=false,
) where {F} # Needed for specialization for this function (avoids some allocations)
    if parallel
        output = map_pairwise_parallel!(
            f,output,box,cl;
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
    map_pairwise!(f::Function,output,box::Box,cl::CellListPair)

The same but to evaluate some function between pairs of the particles of the vectors.

"""
function map_pairwise!(f::F1, output, box::Box, cl::CellListPair{V,N,T,Swap};
    # Parallelization options
    parallel::Bool=true,
    output_threaded=nothing,
    reduce::F2=reduce,
    show_progress::Bool=false
) where {F1,F2,V,N,T,Swap} # F1, F2 Needed for specialization for these functions
    if Swap == Swapped
        fswap(x,y,i,j,d2,output) = f(y,x,j,i,d2,output) 
    else
        fswap = f
    end
    if parallel
        output = map_pairwise_parallel!(
            fswap,output,box,cl;
            output_threaded=output_threaded,
            reduce=reduce,
            show_progress=show_progress
        )
    else
        output = map_pairwise_serial!(fswap,output,box,cl,show_progress=show_progress)
    end
    return output
end

"""
    map_pairwise(args...;kargs...) = map_pairwise!(args...;kargs...)

is an alias for `map_pairwise!` which is defined for two reasons: first, if the output of the funciton is immutable, it may be 
clearer to call this version, from a coding perspective. Second, the python interface through `juliacall` does not accept the 
bang as a valid character. 

"""
const map_pairwise = map_pairwise!

# Utils
include("./neighborlists.jl")

#
# Test and example functions
#
include("./examples/examples.jl")
include("./testing.jl")

# Precompilation tools
include("precompile.jl")

end # module



