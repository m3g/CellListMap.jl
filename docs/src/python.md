# Calling from Python

Callling `CellListMap` from python can be useful if lists of neighbors or other properties have to be computed many times, making the overhead of initializing Julia negligible. As the example and benchmark below demonstrates, the current implementation of cell lists in this package is faster than common alternatives available in the python ecosystem. 

## Installing

First, install `juliacall` using the `pip` package manager, with
```bash
% pip install juliacall
```

Using `ipython3` (only Python $\geq$ 3 is supported), do:

```python
In [1]: from juliacall import Main as jl
```
which, *on the first use only*, will install the latest stable version of Julia. 

Then, install `CellListMap`, with:

```python
In [2]: jl.Pkg.add("CellListMap")
```

## A Python module

The [CellListMap.py](https://github.com/m3g/CellListMap.jl/blob/main/src/examples/CellListMap.py) 
provides a complete small python module that interfaces the `neighborlist` function of `CellListMap` 
with python, returning `numpy` arrays of indices and distances: 

By saving the file above in a `CellListMap.py` file, within python we just need to do:

```python
In [1]: import CellListMap as cl

In [2]: import numpy as np

In [3]: coords = np.random.random((50_000,3))

In [4]: i_inds, j_inds, d = cl.neighborlist(coords, 0.05)
```

The output `i_inds`, `j_inds` and `d` variables are `numpy` arrays with the indexes of the particles and their distances.

For periodic systems, the `unitcell` must be provided, as uni-dimensional `np.array` (for orthorhombic systems) or a `np.matrix` (for
general periodic boundary conditions). For example: 
```python
In [5]: i_inds, j_inds, d = cl.neighborlist(coords, 0.05, unitcell=np.array([1, 1, 1]))

In [6]: i_inds, j_inds, d = cl.neighborlist(coords, 0.05, unitcell=np.matrix('1 0 0; 0 1 0; 0 0 1'))
```

The `neighborlist_cross` function provided above has a similar syntax, but to compute the neighboring particles of two
independent sets:
```python
In [7]: x = np.random.random((50_000,3))

In [8]: y = np.random.random((50_000,3))

In [9]: i_inds, j_inds, d = cl.neighborlist_cross(x, y, 0.05, unitcell=np.array([1, 1, 1]))
```

!!! note
    The indexes of the particles the `i_inds` and `j_inds` arrays are 0-based, to conform the `numpy` array standard. 

!!! tip
    To run the code multi-threaded, set the `JULIA_NUM_THREADS` environment variable before launching python:
    ```bash
    % export JULIA_NUM_THREADS=8
    ```

## Under the hood: interfacing with the Julia package

!!! note
    The details of the above module are explained below, for a more in depth understanding of the
    interface between Julia and Python through the [`PythonCall.jl`](https://github.com/cjdoris/PythonCall.jl) library.

    We highly recommend using the `CellListMap.py` module provided above.

The typical input coordinates, in python, are a `numpy` array with shape `(N,dim)` where `N` is the number of particles and `dim` is the dimension
of the space (2 or 3 for `CellListMap`). Here, we generate a set of `50,000` particles in three dimensions:
```python
In [1]: import numpy as np

In [2]: coords = np.random.random((50_000,3))
```

Julia is column-major, and python is row-major, thus if we want to use the functions from `CellListMap` we need to transpose the coordinates:
```python
In [3]: coords_t = coords.transpose()
```

These transposed coordinates can be used in the `CellListMap.neighborlist` function. For example:
```python
In [4]: jl.seval("using CellListMap")

In [6]: neighbor_list = jl.neighborlist(coords_t,0.05)
```
which will return a list of tuples, containing all pairs of coordinates withing the cutoff (remember that the *first* call to a Julia function will always take longer than subsequent calls, because the function is JIT compiled):

```python
In [12]: neighbor_list.shape
Out[12]: (618774,)

In [13]: neighbor_list[1]
Out[13]: (1, 37197, 0.047189685889846615)
```
Note that the third element of the tuple is the distance between the points.

### Converting the list to numpy arrays

The output of `CellListMap.neighborlist` is a Julia `Vector{Tuple{Int,Int,Float64}}` array (or `Float32`, if the coordinates
and cutoff were given in 32-bit precision). This Julia list can be accessed from within python normally:
```python
In [36]: neighbor_list = jl.neighborlist(coords_t, 0.05);

In [37]: neighbor_list[0:2]
Out[37]: 
2-element view(::Vector{Tuple{Int64, Int64, Float64}}, 1:1:2) with eltype Tuple{Int64, Int64, Float64}:
 (1, 6717, 0.020052121336342873)
 (1, 7208, 0.03880915662838867)

In [38]: neighbor_list[0][0]
Out[38]: 1

In [40]: neighbor_list[0][2]
Out[40]: 0.020052121336342873
```

Yet, this list may not be interoperable with many other python packages, particularly with `numpy` standard 
operations. Thus, it may be interesting to convert the list to `numpy`  arrays. This can be done with a simple
helper function, which uses a Julia function to copy the list values to the `numpy` arrays:

```python
jl.seval("""
function copy_to_numpy_arrays(nb_list, i_inds, j_inds, d)
    for i in eachindex(nb_list)
        i_inds[i], j_inds[i], d[i] = nb_list[i]
    end
    return nothing
end
""")
def neighborlist(x, cutoff) :
    x_t = x.transpose()
    nb_list = jl.neighborlist(x_t, cutoff)
    i_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    j_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    d = np.full((len(nb_list),), 0.0, dtype=np.float64)
    jl.copy_to_numpy_arrays(nb_list, i_inds, j_inds, d)
    i_inds -= 1 # make indexes 0-based
    j_inds -= 1 # make indexes 0-based
    return i_inds, j_inds, d
```

Now, the output of the python `neighborlist` contains the `numpy` arrays for the indexes
of the two particles involved in each pair, and their distances:
```python
In [61]: neighborlist(coords,0.05)
Out[61]: 
(array([    0,     0,     0, ..., 49802, 49802, 49885]),
 array([ 6717,  7208,  9303, ..., 11542, 27777, 43853]),
 array([0.02005212, 0.03880916, 0.04543936, ..., 0.04671987, 0.02671908,
        0.02772025]))
```

The overhead of these conversions, array creation and copies is not very large, and the benchmarks
below are still valid considering this auxiliary python function.

## Benchmarking vs. Scipy

To properly benchmark the `neighborlist` function from `CellListMap`, let us first define a simple wrapper that will include the transposition of the coordinates in the time:

```python
In [14]: def neighborlist_simple(x,cutoff):
    ...:     y = x.transpose()
    ...:     nn = jl.CellListMap.neighborlist(y,cutoff)
    ...:     return nn
    ...:

In [15]: %timeit neighborlist_simple(coords,0.05)
61.7 ms ± 707 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

Let us compare this with the performance of a [inrange neighborlist algorithm](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_ball_tree.html) from `scipy`:

```python
In [29]: from scipy.spatial import cKDTree

In [30]: def neighborlist_scipy(x,cutoff) : 
    ...:     kd_tree = cKDTree(x)  
    ...:     pairs = kd_tree.query_pairs(r=0.05)  
    ...:     return pairs 
    ...:

In [31]: %timeit neighborlist_scipy(coords,0.05)
312 ms ± 2.85 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Just to confirm, this is the number of pairs that is being output in this test
```python
In [32]: len(neighborlist_scipy(coords,0.05)) # using Scipy
Out[32]: 618475

In [20]: len(neighborlist_smple(coords,0.05)) # using CellListMap
Out[20]: 618475
```

If we use the `neighborlist` function from [Converting the list to numpy arrays](@ref), the result is similar,
thus copying the output to numpy arrays does not create a large overhead:
```python
In [30]: %timeit neighborlist(coords, 0.05)
67.4 ms ± 4.04 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

### Overhead

The overhead of calling the function through `juliacall`  is small. From within Julia, the timings of a similar execution would be:
```julia-repl
julia> using BenchmarkTools

julia> using CellListMap

julia> x = rand(3,50_000);

julia> @btime CellListMap.neighborlist($x,0.05,parallel=false);
  51.299 ms (17687 allocations: 37.43 MiB)
```

### Multi-threading

These examples were run single-threaded. To run multi-threaded, an environment variable for `Julia` needs to be set. For example,
in `bash`, do:
```bash
% export JULIA_NUM_THREADS=12
```

!!! warning 
    There is a conflict between garbage collectors that may cause segmentation faults in multi-threaded runs 
    (see [this issue](https://github.com/cjdoris/PythonCall.jl/issues/201)). The workaround appears to be to 
    disable the Julia garbage collector during the execution of multi-threaded code. 
    
    Here we provide the necessary syntax as an auxiliary Python function.

Consider the following python file, let us call it `neighborlist.py`, that provides the `neighborlist`
python function with the conversion of the output to `numpy` arrays:

```python
from juliacall import Main as jl
jl.seval("using CellListMap")
import numpy as np
jl.seval("""
function copy_to_numpy_arrays(nb_list, i_inds, j_inds, d)
    for i in eachindex(nb_list)
        i_inds[i], j_inds[i], d[i] = nb_list[i]
    end
    return nothing
end
""")
def neighborlist(x, cutoff) :
    x_t = x.transpose()
    jl.GC.enable(False)
    nb_list = jl.neighborlist(x_t, cutoff)
    jl.GC.enable(True)
    i_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    j_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    d = np.full((len(nb_list),), 0.0, dtype=np.float64)
    jl.copy_to_numpy_arrays(nb_list, i_inds, j_inds, d)
    return i_inds, j_inds, d
```

Then, in Python, do:

```python
In [1]: import neighborlist as nb

In [2]: import numpy as np

In [3]: coords = np.random.random((50_000,3))

In [4]: i_inds, j_inds, d = nb.neighborlist(coords, 0.05)
```

In a notebook with 6 cores (12 threads) this led to the following performance:
```python
In [5]: %timeit i_inds, j_inds, d = nb.neighborlist(coords, 0.05)
23.7 ms ± 910 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

Which, is about 3x faster than the serial execution:
```python
In [4]: %timeit i_inds, j_inds, d = nb.neighborlist(coords, 0.05)
59.2 ms ± 959 µs per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

and thus about 10x faster than `scipy.spatial`:
```python
In [7]: %timeit neighborlist_scipy(coords,0.05)
204 ms ± 2.86 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

## General mappings

A greater flexibility on the use of `CellListMap` from python can be obtained by defining custom Julia functions.
This feature must be used with the low level interface of `CellListMap`, and is somewhat limited in scope.

```python
In [36]: jl.seval("using CellListMap")

In [37]: x = np.random.random((50_000,3));

In [38]: x_t = x.transpose()

In [39]: box = jl.Box(np.array([1,1,1]), 0.05)

In [40]: box
Out[41]: 
Box{OrthorhombicCell, 3}
  unit cell matrix = [ 1.0, 0.0, 0.0; 0.0, 1.0, 0.0; 0.0, 0.0, 1.0 ]
  cutoff = 0.05
  number of computing cells on each dimension = [22, 22, 22]
  computing cell sizes = [0.05, 0.05, 0.05] (lcell: 1)
  Total number of cells = 10648

In [41]: cl = jl.CellList(x_t,box)

In [42]: cl
Out[42]: 
CellList{3, Float64}
  50000 real particles.
  7985 cells with real particles.
  66594 particles in computing box, including images.
```

The function to be mapped, however, has to be defined in Julia, using `seval`. For example, here we define a function that computes the histogram of the distances within the cutoff. 

```python
In [43]: jl.seval("""  
    ...: function histogram(pair, hist) 
    ...:     cutoff = 0.05 
    ...:     dc = pair.d/cutoff # in [0,1] 
    ...:     ibin = floor(Int,dc*10) + 1 # in [0,10] 
    ...:     hist[ibin] += 1 
    ...:     return hist 
    ...: end 
    ...: """)
Out[44]: histogram (generic function with 1 method)
```

We can initialize the output variable (the histogram) using a regular `numpy` array: 

```python
In [45]: hist = np.zeros(10)
```

and call the `pairwise` function to obtain the histogram of the distances within the `cutoff`:

```python
In [46]: jl.pairwise!(jl.histogram, hist, box, cl)
Out[46]: 
10-element PythonCall.PyArray{Float64, 1, true, true, Float64}:
 153344.0
      1.151744e6
      3.066624e6
      5.787392e6
      9.220608e6
      1.3175552e7
      1.7414912e7
      2.1817088e7
      2.6189312e7
      3.0583808e7
```

Note that the function receives a `NeighborPair` object (`pair`) with fields `pair.i`, `pair.j`, `pair.x`, `pair.y`, `pair.d2`, and a lazily-computed `pair.d` for the distance.

With this interface, however, it is not possible to pass additional parameters to the mapped function, and thus the additional parameters have to defined inside the called function (as the `cutoff` in the current example). This is not ideal, for example, for computing accelerations, which depend on the masses of the particles. In this case, currently, either just use Julia from start and closures, or use the `neighborlist`  function to obtain the list of neighbors to then compute whatever property is desired from the list of pairs, although this is suboptimal in terms of performance.  




