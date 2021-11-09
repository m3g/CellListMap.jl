# Periodic boundary conditions

## Orthorhombic periodic boundary conditions

Orthorhombic periodic boundary conditions allow some special methods that are faster than those for general cells. To initialize an Orthorhombic cell, just provide the length of the cell on each side, and the cutoff. For example:

```julia-repl
julia> box = Box([100,70,130],12)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [100.0 0.0 0.0; 0.0 70.0 0.0; 0.0 0.0 130.0]
  cutoff: 12.0
  number of computing cells on each dimension: [10, 7, 12]
  computing cell sizes: [12.5, 14.0, 13.0] (lcell: 1)
  Total number of cells: 840
```

## Triclinic periodic boundary conditions

Triclinic periodic boundary conditions of any kind can be used. However, the input has some limitations for the moment. The lattice vectors must have strictly positive coordinates, and the smallest distance within the cell cannot be smaller than twice the size of the cutoff. An error will be produced if the cell does not satisfy these conditions. 

Let us illustrate building a two-dimensional cell, for easier visualization. A matrix of column-wise lattice vectors is provided in the construction of the box, and that is all. 

Here, the lattice vectors are `[1,0]` and `[0.5,1]` (and we illustrate with `cutoff=0.1`): 

```julia-repl
julia> box = Box([ 1.0  0.5
                     0  1.0 ], 0.1);

julia> x = 10*rand(SVector{2,Float64},1000);
```
We have created random coordinates for `1000` particles, that are not necessarily wrapped according to the periodic boundary conditions. We can see the coordinates in the minimum image cell with:
```julia-repl
julia> using Plots

julia> CellListMap.draw_computing_cell(x,box)
```

```@raw html
<img src=https://raw.githubusercontent.com/m3g/CellListMap.jl/main/docs/src/assets/lattice.png>
```

The construction of the cell list is, as always, done with:

```julia-repl
julia> cl = CellList(x,box)
CellList{2, Float64}
  109 cells with real particles.
  2041 particles in computing box, including images.

```

Upon construction of the cell lists, the particles are replicated to fill a rectangular box (or orthorhombic box, in three-dimensions), with boundaries that exceed the actual system size. This improves the performance of the pairwise computations by avoiding the necessity of wrapping coordinates on the main loop (this is an implementation detail only). 

In summary, to use arbitrary periodic boundary conditions, just initialize the box with the matrix of lattice vectors. In three dimensions, for example, one could use:

```julia-repl
julia> box = Box([ 50.  0. 00. 
                    0. 30. 30.          
                    0. 00. 50. ],  2.)

julia> x = 100*rand(SVector{3,Float64},10000);

julia> p = [ CellListMap.wrap_to_first(x,box) for x in x ];

julia> using Plots

julia> scatter(Tuple.(p),aspect_ratio=1,framestyle=:box,label=:none)
```
to work with an arbitrary 3D lattice, Which in this case looks like:

```@raw html
<img src=https://raw.githubusercontent.com/m3g/CellListMap.jl/main/docs/src/assets/3Dlattice.png>
```

## Without periodic boundary conditions

To avoid the use of periodic boundary conditions it is enough to define an Orthorhombic box with lengths in each direction that are larger than the limits of the coordinates of the particles plus the cutoff. This can be done automatically with the `limits` function. The box must be constructed with:

```julia-repl
julia> x = [ [100,100,100] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x),12)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [111.99749159163106 0.0 0.0; 0.0 111.99757156637344 0.0; 0.0 0.0 111.99910298572958]
  cutoff: 12.0
  number of computing cells on each dimension: [11, 11, 11]
  computing cell sizes: [12.444165732403452, 12.444174618485938, 12.444344776192175] (lcell: 1)
  Total number of cells: 1331
```

or, for computing the interaction between two disjoint sets of particles, call the `limits` function with two arguments:

```julia-repl
julia> x = [ [100,100,100] .* rand(3) for i in 1:100_000 ];

julia> y = [ [120,180,100] .* rand(3) for i in 1:100_000 ];

julia> box = Box(limits(x,y),12)
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [131.9978650409108 0.0 0.0; 0.0 191.99730748624336 0.0; 0.0 0.0 111.99917288242698]
  cutoff: 12.0
  number of computing cells on each dimension: [12, 17, 11]
  computing cell sizes: [13.19978650409108, 12.799820499082891, 12.444352542491886] (lcell: 1)
  Total number of cells: 2244
```

Note that the unit cell length is, on each direction, the maximum coordinates of all particles plus the cutoff. This will avoid the computation of pairs of periodic images. The algorithms used for computing interactions in Orthorhombic cells will then be used.




