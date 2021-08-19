"""

```
neighbourlist(box, cl; parallel=true)
```

Compute the neighbour list of a single set of particles. Returns a vector of tuples
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

julia> cl = CellList(x,box)
CellList{3, Float64}
  999 cells with real particles.
  1126 particles in computing box, including images.


julia> CellListMap.neighbourlist(box,cl,parallel=false)
15-element Vector{Tuple{Int64, Int64, Float64}}:
 (187, 511, 0.010346860078531755)
 (203, 708, 0.010777737363239403)
 (296, 579, 0.018124283912224655)
 (356, 715, 0.016284309721945608)
 (362, 656, 0.01919577305081326)
 (407, 944, 0.012943980502242866)
 (463, 379, 0.01897013213107807)
 (500, 793, 0.019053137224533643)
 (530, 780, 0.013460883252038484)
 (544, 367, 0.019006016702941237)
 (558, 225, 0.018190653229807584)
 (584, 4, 0.016935844769524398)
 (725, 749, 0.019971874892397875)
 (773, 119, 0.01835233336121765)
 (927, 8, 0.011234110402648743)

```

"""
function neighbourlist(box, cl; parallel=true)

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
    pairs = Tuple{Int,Int,Float64}[]
  
    # Run pairwise computation
    pairs = map_pairwise!(
      (x, y, i, j, d2, pairs) -> push_pair!(i, j, d2, pairs),
      pairs, box, cl,
      reduce=reduce_pairs, parallel=parallel
    )
    return pairs

end
