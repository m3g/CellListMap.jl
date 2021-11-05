"""

```
neighbourlist(box, cl; parallel=true)
```

Compute the neighbour list of a single set or set pairs of particles. Returns a vector of tuples
with all indices of the particles that are within `box.cutoff`, and the distances.  

### Example
```julia-repl
julia> x = [ rand(3) for i in 1:1000 ];

julia> box = Box([1,1,1],0.02) 
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
  cutoff: 0.02
  number of computing cells on each dimension: [52, 52, 52]
  computing cell sizes: [0.02, 0.02, 0.02] (lcell: 1)
  Total number of cells: 140608

julia> cl = CellList(x,box) # single set
CellList{3, Float64}
  999 cells with real particles.
  1126 particles in computing box, including images.


julia> CellListMap.neighbourlist(box,cl,parallel=false)
15-element Vector{Tuple{Int64, Int64, Float64}}:
 (187, 511, 0.010346860078531755)
 (203, 708, 0.010777737363239403)
 (296, 579, 0.018124283912224655)
 ⋮
 (584, 4, 0.016935844769524398)
 (725, 749, 0.019971874892397875)
 (773, 119, 0.01835233336121765)
 (927, 8, 0.011234110402648743)

```

To obtain the neighbour list (within the cutoff) between two sets of 
particles, initialize the cell lists with the two sets: 

```julia-repl
julia> x = [ rand(3) for i in 1:1000 ];

julia> y = [ rand(3) for i in 1:1000 ];

julia> box = Box([1,1,1],0.02);

julia> cl = CellList(x,y,box)

julia> cl = CellList(x,y,box)
CellListMap.CellListPair{Vector{SVector{3, Float64}}, 3, Float64}
   1000 particles in the reference vector.
   997 cells with real particles of target vector.

julia> CellListMap.neighbourlist(box,cl)
35-element Vector{Tuple{Int64, Int64, Float64}}:
 (409, 982, 0.01634641594779082)
 (521, 422, 0.00919026348035512)
 (625, 731, 0.012986301890746663)
 ⋮
 (647, 730, 0.01565763971458105)
 (296, 668, 0.016556686306217868)
 (992, 589, 0.018392993428289553)

```

"""
function neighbourlist(box::Box, cl; parallel=true)

    # Function adds pair to the list
    function push_pair!(i, j, d2, pairs) 
        d = sqrt(d2)
        push!(pairs, (i, j, d))
        return pairs
    end

    # We have to define our own reduce function here (for the parallel version)
    function reduce_pairs(pairs, pairs_threaded)
        pairs = pairs_threaded[1]
        for i in 2:Threads.nthreads()
            append!(pairs, pairs_threaded[i])
        end
        return pairs
    end
  
    # Initialize
    pairs = Tuple{Int,Int,typeof(box.cutoff)}[]
  
    # Run pairwise computation
    pairs = map_pairwise!(
      (x, y, i, j, d2, pairs) -> push_pair!(i, j, d2, pairs),
      pairs, box, cl,
      reduce=reduce_pairs, parallel=parallel
    )
    return pairs

end

"""

```
neighbourlist(x,r;parallel=true)
```

Computes the list of pairs of particles in `x` which are closer to each other than `r`.

### Example
```julia-repl
julia> x = [ rand(3) for i in 1:10_000 ];

julia> CellListMap.neighbourlist(x,0.05)
24848-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 1055, 0.022977369806392412)
 (1, 5086, 0.026650609138167428)
 ⋮
 (9989, 3379, 0.0467653507446483)
 (9989, 5935, 0.02432728985151653)

```

"""
function neighbourlist(x,r;parallel=true)
    box = Box(limits(x),r)
    cl = CellList(x,box,parallel=parallel)
    return neighbourlist(box,cl,parallel=parallel)
end

"""

```
neighbourlist(x,y,r;parallel=true,autoswap=true)
```

Computes the list of pairs of particles of `x` which are closer than `r` to
the particles of `y`. The `autoswap` option will swap `x` and `y` to try to optimize
the cost of the construction of the cell list. 

### Example
```julia-repl
julia> x = [ rand(3) for i in 1:10_000 ];

julia> y = [ rand(3) for i in 1:1_000 ];

julia> CellListMap.neighbourlist(x,y,0.05)
5006-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 269, 0.04770884036497686)
 (25, 892, 0.03850515231540869)
 ⋮
 (9952, 749, 0.048875643578313456)
 (9984, 620, 0.04101242499363183)

```

"""
function neighbourlist(x,y,r;parallel=true,autoswap=true)
    box = Box(limits(x,y),r)
    cl = CellList(x,y,box,parallel=parallel,autoswap=autoswap)
    return neighbourlist(box,cl,parallel=parallel)
end
