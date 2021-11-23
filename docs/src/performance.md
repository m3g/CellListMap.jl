# Performance tunning and additional options

## Optimizing the cell grid

The partition of the space into cells is dependent on a parameter `lcell` which can be passed to `Box`. For example:
```julia
box = Box(x,box,lcell=2)
cl = CellList(x,box)
map_pairwise!(...)
```
This parameter determines how fine is the mesh of cells. There is a trade-off between the number of cells and the number of particles per cell. For low-density systems, greater meshes are better, because each cell will have only a few particles and the computations loop over a smaller number of cells. For dense systems, it is better to run over more cells with less particles per cell. It is a good idea to test different values of `lcell` to check which is the optimal choice for your system. Usually the best value is `lcell=1`, because in `CellListMap` implements a method to avoid spurious computations of distances on top of the cell lists, but for very dense systems, or for very large cutoffs (meaning, for situations in which the number of particles per cell may be very large), a greater `lcell` may provide a better performance. It is unlikely that `lcell > 3` is useful in any practical situation. For molecular systems with normal densities `lcell=1` is likely the optimal choice. The performance can be tested using the progress meter, as explained below.  

As a rough guide, `lcell > 1` is only worthwhile if the number of particles per cell is greater than  `~200-400`.  

## Output progress 

For long-running computations, the user might want to see the progress. A progress meter can be turned on with the `show_progress` option. For example:
```julia
map_pairwise!(f,output,box,cl,show_progress=true)
```
whill print something like:
```julia-repl
Progress:  43%|█████████████████                    | ETA: 0:18:25
```

Thus, besides being useful for following the progress of a long run, it is useful to test different values of `lcell` to tune the performance of the code, by looking at the estimated time to finish (ETA) and killing the execution after a sample run. The default and recommended option for production runs is to use `show_progress=false`, because tracking the progress introduces a small overhead into the computation. 

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
