# Parallelization splitting and reduction

The parallel execution requires the splitting of the computation among threads, obviously. Thus, the output variable must be split and then reduced to avoid concurrency. To control these steps, set manually the `output_threaded` and `reduce` optional input parameters of the `map_pairwise!` function. 

By default, we define:
```julia
output_threaded = [ deepcopy(output) for i in 1:nbatches(cl) ]
```
where `nbatches(cl)` is the number of batches into which the computation will be divided, as defined for the cell list `cl` (this parameter is by default `2*nthreads()`, but it can be tunned for performance, as explained in the **Number of batches** section below), 
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

Using the `length` of the `output_threaded` vector as the measure of how many copies of the array is available is convenient because it will be insensitive in changes in the number of batches that may be set.

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
forces_threaded = [ deepcopy(forces) for i in 1:nbatches(cl) ]
for i in 1:nsteps
    map_pairwise!(f, forces, box, cl, output_threaded=forces_threaded)
    # work with the final forces vector
    ...
    # Reset forces_threaded
    for i in 1:nbatches(cl)
        @. forces_threaded[i] = zero(SVector{3,Float64}) 
    end
end
```
In this case, the `forces` vector will be updated by the default reduction method. `nbatches(cl)` is the number of batches of the parallel calculation, which is defined on the construction of the cell list (by default twice the number of threads available, see the next section).

## Number of batches

Every calculation with cell lists has two steps: the construction of the lists, and the mapping of the computation among the pairs of particles that satisfy the cutoff criterion. 

The construction of the cell list is harder to parallelize, because assigning each particle to a cell is fast, such that the cost of merging a set of lists generated in parallel can be as costly as building the lists themselves. Therefore, it is frequent that it is not worthwhile (actually it is detrimental for performance) to split the construction of the cell lists in too many threads. This is particularly relevant for smaller systems, for which the cost of constructing the lists can be comparable to the cost of actually computing the mapped function. 

At the same time, the homogeneity of the computation of the mapped function may be fast or not, homogeneous or not. These characteristics affect the optimal workload splitting strategy. For very large systems, or systems for which the function to be computed is not homogeneous in time, it may be interesting to split the workload in many tasks as possible, such that slow tasks do not dominate the final computational time.   

Both the above considerations can be used to tunning the `nbatches` parameter of the cell list. This parameter assumes a value of type `NumberOfBatches`, which is basically a tuple of two integers, defining the number of batches that will be used for constructing the cell lists and for the mapping of the computations. 

By default, the number of batches for the computation of the cell lists is smaller than `nthreads()` if the number of particles per cell is small, and cannot be greater than `nthreads()`. The default value by the internal function `CellListMap._nbatches_build_cell_lists(cl::CellList)`. The default value for the number of batches of the function mapping is `2*nthreads()` for computations involving one set of particles, and `length(x÷2500)` for computations involving two sets of particles, where `x` is the set with the greater number of particles (over which the calculation will be split into threads). 

The values assumed for each number of batches can bee seen by printing the `nbatches` parameter of the cell lists:
```julia-repl
julia> Threads.nthreads()
8

julia> x = [ rand(3) for _ in 1:10_000 ]; box = Box([1,1,1],0.1);

julia> cl = CellList(x,box);

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 2
  Number of batches for function mapping: 64
```
which means that the construction of the cell lists will use 2 batches (thus using less tan `nthreads()` tasks), and the mapping of the function will be split into 64 batches. Using more batches than threads for the function mapping is effective most times in avoiding uneven workload, but it may be a problem if the output to be reduced is too large, as the threaded version of the output contains `nbatches` copies of the output. 

The effect of the number of batches in the construction of the cell lists can be seen here (in the above example, with `10_000` particles):

```julia-repl
julia> @btime CellList($x,$box,nbatches=NumberOfBatches(1,64));
  2.231 ms (4189 allocations: 2.12 MiB)

julia> @btime CellList($x,$box,nbatches=NumberOfBatches(2,64)); # default
  1.437 ms (10472 allocations: 3.30 MiB)

julia> @btime CellList($x,$box,nbatches=NumberOfBatches(3,64)); 
  1.403 ms (14089 allocations: 4.09 MiB)

julia> @btime CellList($x,$box,nbatches=NumberOfBatches(4,64)); 
  1.504 ms (18484 allocations: 4.94 MiB)

julia> @btime CellList($x,$box,nbatches=NumberOfBatches(8,64)); 
  2.699 ms (39864 allocations: 8.42 MiB)
```
and, as shown, the default splitting is close to optimal, even if using less then the number of threads available. The optimal number of batches is, however, problem dependent, and the default heuristic may not always choose the best value.


For denser systems the optimal number of batches change. For example, for `1_000_000` particles, we have:
```julia-repl
julia> @btime CellList($x,$box,nbatches=NumberOfBatches(2,64));
  126.990 ms (8980 allocations: 194.95 MiB)

julia> @btime CellList($x,$box,nbatches=NumberOfBatches(8,64)); # default
  78.814 ms (41038 allocations: 353.82 MiB)
```
the default value is again close to optimal and can be trusted.

The number of batches for the mapping of the pairwise computation generally is optimal if greater than the number of threads. 
Most times it doesn't really makes sense to start a number of batches that is not a multiple of the number of threads available. For example, if the number of batches is `nthreads()+1`, most likely `nthreads()` batches will finish almost simultaneously and then the remaining batch will start running. Let us see the effect  of the number of batches in one specific example. The computer in which these tests are performed has 4 physical cores, which with multi-threading can span 8 independent threads.
```julia-repl
julia> Threads.nthreads()
8
```

The test will be the computation of pairwise velocities of a set of `100_000` particles. We will keep the first parameter fixed (the number of batches of the cell list construction):
```julia-repl
julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=NumberOfBatches(2,16)); # default
  42.208 ms (88779 allocations: 36.13 MiB)

julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=NumberOfBatches(2,8)); # = nthreads()
  42.197 ms (88325 allocations: 36.04 MiB)

julia> @btime CellListMap.Examples.pairwise_velocities(N=100_000,nbatches=NumberOfBatches(2,4)); # < nthreads()
  52.379 ms (88291 allocations: 36.04 MiB)
```
As shown above, the optimal number o batches is close to the number of threads available, and increasing it further does not improve performance. It may degrade performance for much larger number of batches, depending on the system size. However, if the computations where heterogeneous and the cost of each batch is much larger than the cost of spawning the threads, splitting into more batches than threads may be worthwhile. Gains in performance for very large systems are expected with this strategy, thus it is the default behavior. 

Finally, the number of batches is set *on the construction of the cell list*, using the `nbatches` keyword parameter. For example:
```julia-repl
julia> cl = CellList(x,box,nbatches=NumberOfBatches(1,4))
CellList{3, Float64}
  1000000 real particles.
  1000 cells with real particles.
  1727449 particles in computing box, including images.

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 1
  Number of batches for function mapping: 4
```
fine tunning of the performance for a specific problem can be obtained by adjusting this parameter. 

If the number of batches is set as zero for any of the two options, the default value is retained. For example:

```julia-repl
julia> cl = CellList(x,box,nbatches=NumberOfBatches(0,4));

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 8
  Number of batches for function mapping: 4

julia> cl = CellList(x,box,nbatches=NumberOfBatches(4,0));

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 4
  Number of batches for function mapping: 64
```