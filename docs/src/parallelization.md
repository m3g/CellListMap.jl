# Parallelization splitting and reduction

The parallel execution requires the splitting of the computation among threads, obviously. Thus, the output variable must be split and then reduced to avoid concurrency. To control these steps, set manually the `output_threaded` and `reduce` optional input parameters of the `map_pairwise!` function. 

By default, we define (here `using Base.Threads`):
```julia
output_threaded = [ deepcopy(output) for i in 1:nthreads() ]
```
and, for scalars and vectors, the reduction is just the sum of the output per thread:
```julia
reduce(output::Number,output_threaded) = sum(output_threaded)
function reduce(output::Vector,output_threaded) 
    for i in 1:nthreads()
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
    for i in 2:Threads.nthreads()
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

The arrays containing the cell lists can be initialized only once, and then updated. This is useful for iterative runs. Note that, since the list size depends on the box size and cutoff, if the box properties changes some arrays might be increased (never shrinked) on this update. 

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
forces_threaded = [ deepcopy(forces) for i in 1:nthreads() ]
for i in 1:nsteps
    map_pairwise!(f, forces, box, cl, output_threaded=forces_threaded)
    # work with the final forces vector
    ...
    # Reset forces_threaded
    for i in 1:Threads.nthreads()
        @. forces_threaded[i] = zero(SVector{3,Float64}) 
    end
end
```
In this case, the `forces` vector will be updated by the default reduction method.