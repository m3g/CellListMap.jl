# Neighbor lists

- [The `neighborlist` function](@ref)
- [In-place computation of neighbor lists](@ref)

## The `neighborlist` function

The package provides a `neighborlist` function that implements this calculation. Without periodic boundary conditions, just do:

```julia-repl
julia> x = [ rand(2) for _ in 1:10_000 ];

julia> neighborlist(x,0.05)
376457-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 363, 0.04855594810064624)
 (1, 513, 0.03356381123125866)
 (1, 1209, 0.005159666709130686)
 ⋮
 (6575, 7378, 0.03791567990447959)
 (7378, 3450, 0.01748757015908321)
```
or `neighborlist(x,y,r)` for computing the lists of pairs of two sets closer than `r`.

The returning array contains tuples with the index of the particle in the first vector, the index of the particle in the second vector, and their distance.

If periodic boundary conditions are used, the `unitcell` can be provided explicitly as keyword parameters:
```julia-repl
julia> x = [ rand(2) for _ in 1:10_000 ]; 

julia> neighborlist(x, 0.05; unitcell=[1,1])
392100-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 5, 0.03445098850037766)
 (1, 393, 0.039448810592487206)
 (1, 1632, 0.02276457565643465)
 ⋮
 (9501, 9781, 0.03351665194098955)
 (9501, 5429, 0.04199258248973222)
```

In the example above, an `Orthorhombic` cell was assumed, and thus a vector of sides was provided. For general
periodic boundary conditions, a unit cell matrix can be provided, for example:

```julia-repl
julia> neighborlist(x, 0.05; unitcell=[1.0 0.5; 0.5 1.0])
580693-element Vector{Tuple{Int64, Int64, Float64}}:
 (1, 457, 0.03935441952786555)
 (1, 1467, 0.033407692174569875)
 (1, 1767, 0.04490555313598093)
 ⋮
 (3652, 8475, 0.04721628783510375)
 (6260, 8475, 0.04946130971686825)
```





