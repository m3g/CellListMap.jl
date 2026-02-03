"""
    resize_output!(sys::AbstractParticleSystem, n::Int)

Resizes the output array and the auxiliary output arrays used
for multithreading, if the number of particles of the system changed.

The function will error if `Base.resize!` is not defined for the 
type of `system.output`. In this case, a `Base.resize!` method
must be implemented by the user. 

!!! warn
    This function *must* be used whenever the output is dependent on
    the number of particles, and that changes, because it adjust the
    size of the copies of the output variable used for multi-threading.

"""
function resize_output!(sys::AbstractParticleSystem, n::Int)
    resize!(sys.output, n)
    for i in eachindex(sys._output_threaded)
        resize!(sys._output_threaded[i], n)
    end
    return sys
end

#
# Function used to update the properties of the systems
#
"""
    update_unitcell!(system, unitcell)

Function to update the unit cell of the system. The `unicell` must be of the 
same type (`OrthorhombicCell`, `TriclinicCell`) of the original `system` 
(changing the type of unit cell requires reconstructing the system).

The `unitcell` can be a `N×N` matrix or a vector of dimension `N`, where
`N` is the dimension of the system (2D or 3D).

This function can be used to update the system geometry in iterative schemes,
where the size of the simulation box changes during the simulation.

!!! note
    Manual updating of the unit cell of non-periodic systems is not allowed.

# Example

```jldoctest ;filter = r" +Parallelization.*" => ""
julia> using CellListMap, StaticArrays, PDBTools

julia> xpositions = coor(read_pdb(CellListMap.argon_pdb_file));

julia> sys = ParticleSystem(
           xpositions = xpositions,
           unitcell=[21,21,21],
           cutoff = 8.0,
           output = 0.0
       );

julia> update_unitcell!(sys, [30.0, 30.0, 30.0])
ParticleSystem1{output} of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 30.0 0.0 0.0; 0.0 30.0 0.0; 0.0 0.0 30.0 ]
      cutoff = 8.0
      number of computing cells on each dimension = [6, 6, 6]
      computing cell sizes = [10.0, 10.0, 10.0] (lcell: 1)
      Total number of cells = 216
    CellListMap.CellList{3, Float64}
      100 real particles.
      8 cells with real particles.
      800 particles in computing box, including images.
    Parallelization auxiliary data set for 4 batch(es).
    Type of output variable (output): Float64

```

"""
function update_unitcell!(sys, unitcell)
    if unitcelltype(sys) == NonPeriodicCell
        throw(
            ArgumentError(
                """\n
                    Manual updating of the unit cell of non-periodic systems is not allowed.

                """
            )
        )
    end
    sys._box = update_box(sys._box; unitcell)
    return sys
end

"""
    update_cutoff!(system, cutoff)

Function to update the `cutoff`` of the system. 

This function can be used to update the system geometry in iterative schemes.

# Example

Here we initialize a particle system with a cutoff of `8.0` and then update
the cutoff to `10.0`. 

```jldoctest ; filter = r"( +Parallelization.*|CellListMap[.])" => ""
julia> using CellListMap, PDBTools

julia> x = coor(read_pdb(CellListMap.argon_pdb_file));

julia> sys = ParticleSystem(
           xpositions = x,
           unitcell=[21.0,21.0,21.0],
           cutoff = 8.0,
           output = 0.0
       );

julia> update_cutoff!(sys, 10.0)
ParticleSystem1{output} of dimension 3, composed of:
    Box{CellListMap.OrthorhombicCell, 3}
      unit cell matrix = [ 21.0 0.0 0.0; 0.0 21.0 0.0; 0.0 0.0 21.0 ]
      cutoff = 10.0
      number of computing cells on each dimension = [5, 5, 5]
      computing cell sizes = [10.5, 10.5, 10.5] (lcell: 1)
      Total number of cells = 125
    CellListMap.CellList{3, Float64}
      100 real particles.
      8 cells with real particles.
      800 particles in computing box, including images.
    Parallelization auxiliary data set for 4 batch(es).
    Type of output variable (output): Float64
```
"""
function update_cutoff!(sys::ParticleSystem1, cutoff)
    if unitcelltype(sys) == NonPeriodicCell
        sys._box = Box(limits(sys.xpositions), cutoff)
    end
    sys._box = update_box(sys._box; cutoff)
    return sys
end
function update_cutoff!(sys::ParticleSystem2, cutoff)
    if unitcelltype(sys) == NonPeriodicCell
        sys._box = Box(limits(sys.xpositions, sys.ypositions), cutoff)
    end
    sys._box = update_box(sys._box; cutoff)
    return sys
end