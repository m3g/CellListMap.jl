# Periodic boundary conditions

Triclinic periodic boundary conditions of any kind can be used. However, the input has some limitations for the moment. The lattice vectors must have strictly positive coordinates, and the smallest distance within the cell cannot be smaller than twice the size of the cutoff. An error will be produced if the cell does not satisfy these conditions. 

Let us illustrate building a two-dimensional cell, for easier visualization. A matrix of column-wise lattice vectors is provided in the construction of the box, and that is all. 

Here, the lattice vectors are `[1,0]` and `[0.5,1]` (and we illustrate with `cutoff=0.1`): 

```julia
julia> box = Box([ 1.0  0.5
                     0  1.0 ], 0.1);

julia> x = 10*rand(SVector{2,Float64},1000);
```
We have created random coordinates for `1000` particles, that are not necessarily wrapped according to the periodic boundary conditions. We can see the coordinates in the minimum image cell with:
```julia
julia> using Plots

julia> CellListMap.draw_computing_cell(x,box)
```

<img src=./src/assets/lattice.png>

The construction of the cell list is, as always, done with:

```julia
julia> cl = CellList(x,box)
CellList{2, Float64}
  109 cells with real particles.
  2041 particles in computing box, including images.

```

Upon construction of the cell lists, the particles are replicated to fill a rectangular box (or orthorhombic box, in three-dimensions), with boundaries that exceed the actual system size. This improves the performance of the pairwise computations by avoding the necessity of wrapping coordinates on the main loop (this is an implementation detail only). 

In summary, to use arbitrary periodic boundary conditions, just initialize the box with the matrix of lattice vectors. In three dimensions, for example, one could use:

```julia
julia> box = Box([ 50.  0. 00. 
                    0. 30. 30.          
                    0. 00. 50. ],  2.)

julia> x = 100*rand(SVector{3,Float64},10000);

julia> p = [ CellListMap.wrap_to_first(x,box) for x in x ];

julia> scatter(Tuple.(p),aspect_ratio=1,framestyle=:box,label=:none)
```
to work with an arbitrary 3D lattice, Which in this case looks like:

<img src=./src/assets/3Dlattice.png>
