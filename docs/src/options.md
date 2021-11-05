# Additional options

## Input coordinates as matrices

For compatibility with other software, the input coordinates can be provided as matrices. The matrices must have dimensions `(2,N)` or `(3,N)`, where `N` is the number of particles (because Julia is column-major, thus this has the same memory layout of an array of length `N` of static vectors). 

For example:
```julia-repl
julia> x = rand(3,100);

julia> box = Box([1,1,1],0.1);

julia> cl = CellList(x,box)
CellList{3, Float64}
  100 real particles.
  99 cells with real particles.
  162 particles in computing box, including images.

julia> map_pairwise!((x,y,i,j,d2,n) -> n += 1, 0, box, cl) # count neighbours
23
```







