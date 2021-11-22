# Parallelization splitting and reduction

The parallel execution requires the splitting of the computation among threads, obviously. Thus, the output variable must be split and then reduced to avoid concurrency. To control these steps, set manually the `output_threaded` and `reduce` optional input parameters of the `map_pairwise!` function. 

By default, we define:
```julia
output_threaded = [ deepcopy(output) for i in 1:CellListMap.batches(cl) ]
```
where `CellListMap.batches(cl)` is the number of batches into which the computation will be divided, as defined for the cell list `cl` (this parameter is usually equal to the number of threads available, except for very small system, but it can be tunned for performance, as explained [here](voltar)), 
and, for scalars and vectors, the reduction is just the sum of the output per thread:
```julia
reduce(output::Number,output_threaded) = sum(output_threaded)
function reduce(output::Vector,output_threaded) 
    @. output = output_threaded[1]
    for i in 2:length(output_threaded)
         @. output += output_threaded[i] 
    end
    return output
end
```

## Custom reduction functions

In some cases, as in the [Nearest neighbour](#nearest-neighbour) example, the output is a tuple and reduction consists in keeping the output from each thread having the minimum value for the distance. Thus, the reduction operation is not a simple sum over the elements of each threaded output. We can, therefore, overwrite the default reduction method, by passing the reduction function as the `reduce` parameter of `map_pairwise!`:
```julia
mind = map_pairwise!( 
    (x,y,i,j,d2,mind) -> f(i,j,d2,mind), mind,box,cl;
    reduce=reduce_mind
)
```
where here the `reduce` function is set to be the custom function that keeps the tuple associated to the minimum distance obtained between threads:
```julia
function reduce_mind(output,output_threaded)
    mind = output_threaded[1]
    for i in 2:length(output_threaded)
        if output_threaded[i][3] < mind[3]
            mind = output_threaded[i]
        end
    end
    return mind
end
```
This function *must* return the updated `output` variable, being it mutable or not, to be compatible with the interface.  

## Preallocating auxiliary arrays: threaded output and cell lists

### Preallocating the cell lists and cell list auxiliary arrays

The arrays containing the cell lists can be initialized only once, and then updated. This is useful for iterative runs. Note that, since the list size depends on the box size and cutoff, if the box properties changes some arrays might be increased (never shrink) on this update. 

```julia
# Initialize cell lists with initial coordinates
cl = CellList(x,box)
# Allocate auxiliary arrays for threaded cell list construction
aux = CellListMap.AuxThreaded(cl)
for i in 1:nsteps
    x = ... # new coordinates
    box = Box(sides,cutoff) # perhaps the box has changed
    cl = UpdateCellList!(x,box,cl,aux) 
end
```

The procedure is identical if using two sets of coordinates, in which case, one would do:
```julia
cl = CellList(x,y,box)
aux = CellListMap.AuxThreaded(cl)
for i in 1:nsteps
    x = ... # new coordinates
    box = Box(sides,cutoff) # perhaps the box has changed
    cl = UpdateCellList!(x,y,box,cl,aux)
end
```

By passing the `aux` auxiliary structure, the `UpdateCellList!` functions will only allocate some minor variables associated to the launching of multiple threads and, possibly, to the expansion of the cell lists if the box or the number of particles became greater. 

### Preallocating threaded output auxiliary arrays

On parallel runs, note that `output_threaded` is, by default, initialized on the call to `map_pairwise!`. Thus, if the calculation must be run multiple times (for example, for several steps of a trajectory), it is probably a good idea to preallocate the threaded output, particularly if it is a large array. For example, the arrays of forces should be created only once, and reset to zero after each use:
```julia
forces = zeros(SVector{3,Float64},N)
forces_threaded = [ deepcopy(forces) for i in 1:cl.nbatches ]
for i in 1:nsteps
    map_pairwise!(f, forces, box, cl, output_threaded=forces_threaded)
    # work with the final forces vector
    ...
    # Reset forces_threaded
    for i in 1:cl.nbatches
        @. forces_threaded[i] = zero(SVector{3,Float64}) 
    end
end
```
In this case, the `forces` vector will be updated by the default reduction method. `cl.nbatches` is the number of batches of the parallel calculation, which is defined on the construction of the cell list (usually equal to the number of threads available).

## Number of batches

Every calculation with cell lists has two steps: the construction of the lists, and the mapping of the computation among the pairs of particles that satisfy the cutoff criterion. 

The construction of the cell list is harder to parallelize, because assigning each particle to a cell is fast, such that the cost of merging a set of lists generated in parallel can be as costly as building the lists themselves. Therefore, it is frequent that it is not worthwhile (actually it is detrimental for performance) to split the construction of the cell lists in to many threads. This is particularly relevant for smaller systems, for which the cost of constructing the lists can be comparable to the cost of actually computing the mapped function. 

At the same time, the homogeneity of the computation of the mapped function may be fast or not, homogeneous or not. These characteristics affect the optimal workload splitting strategy. For very large systems, or systems for which the function to be computed is not homogeneous in time, it may be interesting to split the workload in many tasks as possible, such that slow tasks do not dominate the final computational time.   

Both the above considerations can be taken into consideration by tunning the `nbatches` parameter in the construction of the cell lists. This parameter assumes a value of type `NumberOfBatches`, which is basically a tuple of two integers, defining the number of batches that will be used for constructing the cell lists and for the mapping of the computations. By default, these assume values which are at most `nthreads()`, but in particular for very small systems the number of batches for the construction of the cell lists can be smaller.

For example:
```julia-repl
julia> Threads.nthreads()
8

julia> x = [ rand(3) for _ in 1:1000 ]; box = Box([1,1,1],0.1);

julia> cl = CellList(x,box); # default

julia> cl.nbatches
CellListMap.NumberOfBatches
  Number of batches for cell list construction: 1
  Number of batches for function mapping: 8
```

Note that we have 8 threads  available, but by default, for this small system (1000 particles), we will use only one batch for the construction of the cell lists. 
This is effectively faster than if we force the computation to use all threads:
```
julia> @btime CellList($x,$box); # default
  262.040 μs (2224 allocations: 627.38 KiB)

julia> @btime CellList($x,$box,nbatches=CellListMap.NumberOfBatches(8,8));
  2.740 ms (18883 allocations: 2.86 MiB)
```

For larger systems this may change. For example, for `1_000_000` points, we have:
```
julia> @btime CellList($x,$box,nbatches=CellListMap.NumberOfBatches(1,8));
  207.017 ms (5227 allocations: 77.14 MiB)

julia> @btime CellList($x,$box,nbatches=CellListMap.NumberOfBatches(8,8));
  78.663 ms (40540 allocations: 351.79 MiB)
```

The default numbers of batches, in these cases, correspond to the optimal choices. But what is implemented is a very simple heuristic based on the number of particles per cell, thus the optimal choice is available for exploration by the user by adjusting the `nbatches` parameter as shown above. 

The number of batches for the mapping of the pairwise computation generally is close to optimal if equal to the number of threads. 
Most times it doesn't really makes sense to start a number of batches that is not a multiple of the number of threads available. For example, if the number of batches is `nthreads()+1`, most likely `ntreads()` batches will finish almost simultaneously and then the remaining batch will start running. Let us see the effect  of the number of batches in one specific example. The computer in which these tests are performed has 4 physical cores, which with multi-threading can span 8 independent threads.
```julia-repl
julia> Threads.nthreads()
8
```

The test will be the computation of pairwise velocities of a set of `100_000` particles. We will keep the first parameter fixed (the number of batches of the cell list construction):
```julia-repl
julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=CellListMap.NumberOfBatches(4,1));
  159.174 ms (178841 allocations: 51.94 MiB)

julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=CellListMap.NumberOfBatches(4,4));
  56.735 ms (178872 allocations: 51.95 MiB)

julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=CellListMap.NumberOfBatches(4,8));
  46.302 ms (178914 allocations: 51.95 MiB)

julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=CellListMap.NumberOfBatches(4,16));
  46.211 ms (178998 allocations: 51.97 MiB)

julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000);
  48.192 ms (88325 allocations: 36.04 MiB)

```
As shown above, the optimal number o batches is close to the number of threads available, and increasing it further does not improve performance. It may degrade performance for larger number of batches. However, if the computations where heterogeneous and the cost of each batch is much larger than the cost of spawning the threads, splitting into more batches than threads may be worthwhile. Gains in performance for very large systems are expected with this strategy. 

For example, for a system with 5 million particles, this computation is faster with more batches than with the number of threads:
```
julia> @time CellListMap.Examples.pairwise_velocities(N=5_000_000,nbatches=CellListMap.NumberOfBatches(4,16));
  3.707087 seconds (8.36 M allocations: 2.066 GiB, 4.86% gc time)

julia> @time CellListMap.Examples.pairwise_velocities(N=5_000_000,nbatches=CellListMap.NumberOfBatches(4,8));
  4.450609 seconds (8.36 M allocations: 2.066 GiB, 20.86% gc time)
```
this occurs because the random splitting of the workload into batches can lead to fluctuations in the amount of work per batch, and this effect is minimized if a greater number of batches is present. Also, this alleviates the problem of one batch being assigned to a thread that for some reason became slower (for competing with other process, or being associated to a virtual processor, etc.)


