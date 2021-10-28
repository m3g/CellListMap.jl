# Performance tunning and additional options

## Optimizing the cell grid

Until this is automatized (hopefuly soon), the partition of the space into cells is dependent on a parameter `lcell` which can be passed to `Box`. For example:
```julia
box = Box(x,box,lcell=2)
cl = CellList(x,box)
map_pairwise!(...)
```
This parameter determines how fine is the mesh of cells. There is a trade-off between the number of cells and the number of particles per cell. For low-density systems, greater meshes are better, because each cell will have only a few particles and the computations loop over a samller number of cells. For dense systems, it is better to run over more cells with less particles per cell. It is a good idea to test different values of `lcell` to check which is the optimal choice for your system. Usually the best value is between `lcell=1` and `lcell=6`, but for large and dense systems a larger value may be optimal. For molecular systems with normal densities `lcell=1` is likely the optimal choice. The peformance can be tested using the progress meter, as explained below.  

## Output progress 

For long-running computations, the user might want to see the progress. A progress meter can be turned on with the `show_progress` option. For example:
```julia
map_pairwise!(f,output,box,cl,show_progress=true)
```
whill print something like:
```julia-repl
Progress:  43%|█████████████████                    | ETA: 0:18:25
```

Thus, besides being useful for following the progress of a long run, it is useful to test different values of `lcell` to tune the peformance of the code, by looking at the estimated time to finish (ETA) and killing the execution after a sample run. The default and recommended option for production runs is to use `show_progress=false`, because tracking the progress introduces a small overhead into the computation. 

## Some benchmarks

### Computing a histogram of pairwise velocities

The goal here is to provide a good implementation of cell lists. We compare it with the implementation of the nice cython/python [halotools](https://github.com/astropy/halotools) package, in the computation of an histogram of mean pairwise velocities. 

```@raw html
<center>
<img src=https://raw.githubusercontent.com/lmiq/PairVelocities/main/data/cd_v0.5.3.png>
<br>
<img src=https://raw.githubusercontent.com/lmiq/PairVelocities/main/data/cv_v0.5.3.png>
</center>
```

The full test is available [at this](https://github.com/lmiq/PairVelocities) repository. And we kindly thank [Carolina Cuesta](https://github.com/florpi) for providing the example. These benchmarks were run on an Intel i7 8th gen laptop, with 4 cores (8 threads). 

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

## Non-allocating type conversion 

Internally, `CellListMap` works with static arrays  for optimal performance. However, sometimes the coordinates are provided using some custom type for which conversion to static arrays is not defined. For example, it is common to use coordinates with units:

```julia-repl
julia> using Unitful

julia> x = [ rand(3)u"nm" for i in 1:100 ];

julia> x[1]
3-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 0.5555681432511039 nm
 0.3112134334494392 nm
 0.6849761663523335 nm

 ```

 In order to use the type of coordinates without allocations and complications in `CellListMap`, just overload the `CellListMap.strip_value` function such that it converts a value of the given type to a float. For example, the `Unitful` package provides the `ustrip` function for that. We define, then:

 ```julia-repl
 julia> CellListMap.strip_value(x::Quantity) = Unitful.ustrip(x)

 ```

 such that it converts a single value of type `Quantity` to a standard float:
```julia-repl
julia> CellListMap.strip_value(x[1][1])
0.5555681432511039

```

With that, the `Unitful` quantities can be passed to `CellList` without modification:

```julia-repl
julia> box = Box([1,1,1],0.1);

julia> CellList(x,box)
CellList{3, Float64}
  100 real particles.
  93 cells with real particles.
  167 particles in computing box, including images.

```








