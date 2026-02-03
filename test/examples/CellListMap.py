#
# This is a Python module providing an interface for CellListMap, particularly
# to the neighborlist functions
# 
# see: https://m3g.github.io/CellListMap.jl/stable/python/
#
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
# Neighboring pairs of a single set of particles
def neighborlist(x, cutoff, unitcell=None) :
    x_t = x.transpose()
    # Garbage collection must be turned off because it does not 
    # interact nicely with python objects, which may be cleaned up
    # causing memory corruption
    jl.GC.enable(False)
    nb_list = jl.neighborlist(x_t, cutoff, unitcell=unitcell)
    jl.GC.enable(True)
    i_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    j_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    d = np.full((len(nb_list),), 0.0, dtype=np.float64)
    jl.copy_to_numpy_arrays(nb_list, i_inds, j_inds, d)
    # convert to 0-based indexing
    i_inds -= 1
    j_inds -= 1
    # Run julia garbage collection
    jl.seval("nb_list = nothing")
    jl.GC.gc()
    return i_inds, j_inds, d
# Neighboring pairs of two sets of particles
def neighborlist_cross(x, y, cutoff, unitcell=None) :
    x_t = x.transpose()
    y_t = y.transpose()
    # Garbage collection must be turned off because it does not 
    # interact nicely with python objects, which may be cleaned up
    # causing memory corruption
    jl.GC.enable(False)
    nb_list = jl.neighborlist(x_t, y_t, cutoff, unitcell=unitcell)
    jl.GC.enable(True)
    i_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    j_inds = np.full((len(nb_list),), 0, dtype=np.int64)
    d = np.full((len(nb_list),), 0.0, dtype=np.float64)
    jl.copy_to_numpy_arrays(nb_list, i_inds, j_inds, d)
    # convert to 0-based indexing
    i_inds -= 1
    j_inds -= 1
    # Run julia garbage collection
    jl.seval("nb_list = nothing")
    jl.GC.gc()
    return i_inds, j_inds, d