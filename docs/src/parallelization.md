# Parallelization splitting and reduction

The parallel execution requires the splitting of the computation among tasks. 

### How output is updated thread-safely

To allow general output types, the approach of `CellListMap` is to copy the output variable the number of times necessary for each parallel task to update an independent output variables, which are reduced at the end. This, of course, requires some additional memory, particularly if the output being updated is formed by arrays. These copies can be preallocated, and custom reduction functions can be defined. 

To control these steps, set manually the `output_threaded` and `reduce` optional input parameters of the `map_pairwise!` function. 

By default, we define:
```julia
output_threaded = [ deepcopy(output) for i in 1:nbatches(cl) ]
```
where `nbatches(cl)` is the number of batches into which the computation will be divided. The number of batches is *not* necessarily equal to the number of threads available (an heuristic is used to optimize performance, as a function of the workload per batch), but can be manually set, as described in the **Number of batches** section below. 

The default reduction function just assumes the additivity of the results obtained by each batch:
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

In some cases, as in the [Nearest neighbor](#nearest-neighbor) example, the output is a tuple and reduction consists in keeping the output from each thread having the minimum value for the distance. Thus, the reduction operation is not a simple sum over the elements of each threaded output. We can, therefore, overwrite the default reduction method, by passing the reduction function as the `reduce` parameter of `map_pairwise!`:
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

## [Number of batches](@id Number-of-batches)

Every calculation with cell lists has two steps: the construction of the lists, and the mapping of the computation among the pairs of particles that satisfy the cutoff criterion. 

The construction of the cell list is harder to parallelize, because assigning each particle to a cell is fast, such that the cost of merging a set of lists generated in parallel can be as costly as building the lists themselves. Therefore, it is frequent that it is not worthwhile (actually it is detrimental for performance) to split the construction of the cell lists in too many threads. This is particularly relevant for smaller systems, for which the cost of constructing the lists can be comparable to the cost of actually computing the mapped function. 

At the same time, the homogeneity of the computation of the mapped function may be fast or not, homogeneous or not. These characteristics affect the optimal workload splitting strategy. For very large systems, or systems for which the function to be computed is not homogeneous in time, it may be interesting to split the workload in many tasks as possible, such that slow tasks do not dominate the final computational time.   

Both the above considerations can be used to tunning the `nbatches` parameter of the cell list. This parameter is initialized from a tuple of integers, defining the number of batches that will be used for constructing the cell lists and for the mapping of the computations. 

By default, the number of batches for the computation of the cell lists is smaller than `nthreads()` if the number of particles per cell is small. The default value by the internal function `CellListMap._nbatches_build_cell_lists(cl::CellList)`. 

The values assumed for each number of batches can bee seen by printing the `nbatches` parameter of the cell lists:
```julia-repl
julia> Threads.nthreads()
64

julia> x, box = CellListMap.xatomic(10^4) # random set with atomic density of water

julia> cl = CellList(x,box);

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 8 
  Number of batches for function mapping: 32 
```
The construction of the cell lists is performed by creating copies of the data, and currently does not scale very well. Thus, no more than 8 batches are used by default, to avoid delays associated to data copying and gargabe collection. The number of batches of the mapping function uses an heuristic which currently limits somewhat the number of batches for small systems, when the overhead of spawning tasks is greater than the computation. 
Using more batches than threads for the function mapping is effective most times in avoiding uneven workload, but it may be a problem if the output to be reduced is too large, as the threaded version of the output contains `nbatches` copies of the output. 

Using less batches than the number of threads also allows the efficient use of nested multi-threading, as the computations will only use the number of threads required, leaving the other threads available for other tasks.

The number of batches is set *on the construction of the cell list*, using the `nbatches` keyword parameter. For example:
```julia-repl
julia> cl = CellList(x,box,nbatches=(1,4))
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
julia> cl = CellList(x,box,nbatches=(0,4));

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 8 
  Number of batches for function mapping: 4

julia> cl = CellList(x,box,nbatches=(4,0));

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 4
  Number of batches for function mapping: 64
```

The number of batches can also be retrieved from the cell list using the `nbatches` function:
```julia-repl
julia> cl = CellList(x,box,nbatches=(2,4));

julia> cl.nbatches
NumberOfBatches
  Number of batches for cell list construction: 2
  Number of batches for function mapping: 4

julia> nbatches(cl) # returns cl.nbatches.map_computation
4

julia> nbatches(cl,:map) # returns cl.nbatches.map_computation
4

julia> nbatches(cl,:build) # returns cl.nbatches.build_cell_lists
2
```

The call `nbatches(cl)` is important for defining the number of copies of preallocated threaded output variables, as explained in the previous section.
