# Neighbor lists

The package provides a `neighborlist` function that implements this calculation. Without periodic boundary conditions, just do:

```julia-repl
julia> x = [ rand(2) for _ in 1:10_000 ];

julia> CellListMap.neighborlist(x,0.05)
24777-element Vector{Tuple{Int64, Int64, Float64}}:
 (0, 62, 0.028481068525796384)
 ⋮
 (9953, 1749, 0.04887502372299809)
 (9973, 124, 0.040110356034451795)
```
or `CellListMap.neighborlist(x,y,r)` for computing the lists of pairs of two sets closer than `r`.

The returning array contains tuples with the index of the particle in the first vector, the index of the particle in the second vector, and their distance.

If periodic boundary conditions are used, the `Box` and `CellList` must be constructed in advance:
```julia-repl
julia> x = [ rand(2) for _ in 1:10_000 ]; 

julia> box = Box([1,1,1],0.1);

julia> cl = CellList(x,box);

julia> CellListMap.neighborlist(box,cl)

julia> CellListMap.neighborlist(box,cl)
209505-element Vector{Tuple{Int64, Int64, Float64}}:
 (0, 121, 0.05553035041478053)
 (0, 1589, 0.051415489701932444)
 ⋮
 (7468, 7946, 0.09760096646331885)
```
