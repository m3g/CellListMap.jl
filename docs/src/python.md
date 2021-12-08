# Calling from Python

Callling `CellListMap` from python can be useful if lists of neighbors or other properties have to be computed many times, making the overhead of initializing Julia negligible. As the example and benchmark below demonstrates, the current implementation of cell lists in this package is faster than common alternatives available in the python ecosystem. 

## Installing

First, install `julia` through the `pip` package manager, with
```bash
% pip install julia
```

Using `ipython3` (only Python $\geq$ 3 is supported), do:

```python
In [1]: from juliacall import Main as jl
```
which, *on the first call only*, will install the latest stable version of Julia. Then, install `CellListMap`, with:
```python
In [2]: jl.Pkg.add("CellListMap")
```

## Calling `neighborlist` 

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
In [4]: from juliacall import Main as jl

In [5]: jl.seval("using CellListMap")

In [6]: neighbor_list = jl.CellListMap.neighborlist(coords_t,0.05)
```
which will return a list of tuples, containing all pairs of coordinates withing the cutoff (remember that the *first* call to a Julia function will always take longer than subsequent calls, because the function is JIT compiled):

```python
In [12]: neighbor_list.shape
Out[12]: (618774,)

In [13]: neighbor_list[1]
Out[13]: (1, 37197, 0.047189685889846615)
```
Note that the third element of the tuple is the distance between the points.

### Benchmark

To properly benchmark the `neighborlist` function from `CellListMap`, let us first define a simple wrapper that will include the transposition of the coordinates in the time:

```python
In [14]: def neighborlist(x,cutoff):
    ...:     y = x.transpose()
    ...:     nn = jl.CellListMap.neighborlist(y,cutoff)
    ...:     return nn
    ...:

In [15]: %timeit neighborlist(coords,0.05)
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
206 ms ± 2.07 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

Just to confirm, this is the number of pairs that is being output in this test
```python
In [32]: len(neighborlist_scipy(x,cutoff)) # using Scipy
Out[32]: 618475

In [20]: len(neighborlist(coords,0.05)) # using CellListMap
Out[20]: 618475
```

## Multi-threading

These examples were run single-threaded. To run multi-threaded, an environment variable for `Julia` needs to be set. For example,
in `bash`, do:
```bash
% export JULIA_NUM_THREADS=8
```

For the current example, this provides a small additional speedup:
```python
In [11]: %timeit neighborlist(coords,0.05)
45.2 ms ± 2.67 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

### Overhead

The overhead of calling the function through `juliacall`  is small. From within Julia, the timings of a similar execution would be:
```julia-repl
julia> using BenchmarkTools

julia> using CellListMap

julia> x = rand(3,50_000);

julia> @btime CellListMap.neighborlist($x,0.05);
  32.786 ms (187997 allocations: 91.03 MiB)

julia> @btime CellListMap.neighborlist($x,0.05,parallel=false);
  51.328 ms (17543 allocations: 32.83 MiB)
```

## General mappings

A greater flexibility on the use of `CellListMap` from python can be obtained by defining custom Julia functions. The construction of the systems and of the cell lists can be performed without modification. For example:

```python
In [28]: box = jl.Box(np.array([1,1,1]),0.05)

In [29]: box

Out[29]: 
Box{OrthorhombicCell, 3, Float64, 9}
  unit cell matrix: [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
  cutoff: 0.05
  number of computing cells on each dimension: [22, 22, 22]
  computing cell sizes: [0.05, 0.05, 0.05] (lcell: 1)
  Total number of cells: 10648
```

```python
In [30]: x = np.random.random((50_000,3))

In [31]: x_t = x.transpose()

In [32]: cl = jl.CellList(x_t,box)

In [33]: cl
Out[33]: 
CellList{3, Float64}
  50000 real particles.
  7982 cells with real particles.
  66532 particles in computing box, including images.
```

The function to be mapped, however, has to be defined in Julia, using `seval`. For example, here we define a function that computes the histogram of the distances within the cutoff. 

```python
In [34]: jl.seval("""  
    ...: function histogram(x,y,i,j,d2,hist) 
    ...:     cutoff = 0.05 
    ...:     dc = sqrt(d2)/cutoff # in [0,1] 
    ...:     ibin = floor(Int,dc*10) + 1 # in [0,10] 
    ...:     hist[ibin] += 1 
    ...:     return hist 
    ...: end 
    ...: """)
Out[34]: histogram (generic function with 1 method)
```

We can initialize the output variable (the histogram) using a regular `numpy` array: 

```python
In [8]: hist = np.zeros(10)
```

and call the `map_pairwise` function to obtain the histogram of the distances within the `cutoff`:

```python
In [37]: jl.map_pairwise(jl.histogram, hist, box, cl)
Out[37]: 
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

With this interface, however, it is not possible to pass additional parameters to the mapped function, and thus the additional parameters have to defined inside the called function (as the `cutoff` in the current example). This is not ideal, for example, for computing accelerations, which depend on the masses of the particles. In this case, currently, either just use Julia from start and closures, or use the `neighborlist`  function to obtain the list of neighbors to then compute whatever property is desired from the list of pairs, although this is suboptimal in terms of performance.  




