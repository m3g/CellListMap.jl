# Types of variables that have support for multi-threading without having
# to explicit add methods to copy_output, reset_output!, and reducer functions.
const SupportedTypes = Union{Number, SVector, FieldVector}

# Supported types for coordinates
const SupportedCoordinatesTypes = Union{Nothing, AbstractVector{<:AbstractVector}, AbstractMatrix}

# Abstract type only for cleaner dispatch
abstract type AbstractParticleSystem{OutputName} end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `pairwise!` computations,
for systems with one set of positions (thus, replacing the loops over `N(N-1)` 
pairs of particles of the set). 

The `xpositions`, `output`, and `parallel` fields are considered part of the API,
and you can retrieve or mutate `xpositions`, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fields of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `ParticleSystem1` structure
is done through the `ParticleSystem(;xpositions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct ParticleSystem1{OutputName, V, O, B, C, A, VC} <: AbstractParticleSystem{OutputName}
    xpositions::V
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
    validate_coordinates::VC
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `pairwise!` computations,
for systems with two set of positions (thus, replacing the loops over `N×M` 
pairs of particles, being `N` and `M` the number of particles of each set).

The `xpositions`, `ypositions`, `output`, and `parallel` fields are considered part of the API,
and you can retrieve or mutate positions, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fields of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `ParticleSystem1` structure
is done through the `ParticleSystem(;xpositions, ypositions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct ParticleSystem2{OutputName, V, O, B, C, A, VC} <: AbstractParticleSystem{OutputName}
    xpositions::V
    ypositions::V
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
    validate_coordinates::VC
end
ParticleSystem2{OutputName}(vx::V, vy::V, o::O, b::B, c::C, vo::Vector{O}, a::A, p::Bool, vc::VC) where {OutputName, V, O, B, C, A, VC} =
    ParticleSystem2{OutputName, V, O, B, C, A, VC}(vx, vy, o, b, c, vo, a, p, vc)

"""
    ParticleSystem(;
        xpositions::Union{AbstractVector{<:AbstractVector},AbstractMatrix},
        #or
        xpositions::Union{AbstractVector{<:AbstractVector},AbstractMatrix},
        ypositions::Union{AbstractVector{<:AbstractVector},AbstractMatrix},
        # and
        unitcell::Union{Nothing,AbstractVecOrMat} = nothing,
        cutoff::Number,
        output::Any;
        output_name::Symbol,
        parallel::Bool=true,
        nbatches::Tuple{Int,Int}=(0, 0),
        validate_coordinates::Union{Nothing,Function}=_validate_coordinates
    )

Constructor of the `ParticleSystem` type given the positions of the particles.

- Positions can be provided as vectors of 2D or 3D vectors 
  (preferentially static vectors from `StaticArrays`), or as 
  (2,N) or (3,N) matrices (v0.8.28 is required for matrices).

- If only the `xpositions` array is provided, a single set of coordinates 
  is considered, and the computation will be mapped for the `N(N-1)` 
  pairs of this set. 

- If the `xpositions` and `ypositions` arrays of coordinates are provided, 
  the computation will be mapped to the `N×M` pairs of particles, 
  being `N` and `M` the number of particles of each set of coordinates.

The unit cell (either a vector for `Orthorhombic` cells or a 
full unit cell matrix for `Triclinic` cells - where columns contain
the lattice vectors), the cutoff used for the
construction of the cell lists and the output variable of the calculations.
If unitcell == nothing, the system is considered not-periodic, in which
case artificial periodic boundaries will be built such that images 
are farther from each other than the cutoff.

!!! note
    The `output` value is the initial value of the output. Tipicaly this is 
    set to `zero(typeof(output))`. In subsequent call to `pairwise!`, 
    the initial value can be optionally reset to `zero(typeof(output))`.

`output_name` can be set to a symbol that best identifies the output variable.
For instance, if `output_name=:forces`, the forces can be retrieved from the
structure using the `system.forces` notation.

The `parallel` and `nbatches` flags control the parallelization scheme of
computations (see https://m3g.github.io/CellListMap.jl/stable/parallelization/#Number-of-batches)).
By default the parallelization is turned on and `nbatches` is set with heuristics
that may provide good efficiency in most cases. 

The `validate_coordinates` function can be used to validate the coordinates
before the construction of the system. If `nothing`, no validation is performed.
By default the validation checks if the coordinates are not missing or NaN. 

# Example

In these examples, we compute the sum of the squared distances between
the particles that are within the cutoff:

## Single set of particles

```jldoctest ;filter = r"(\\d*)\\.(\\d)\\d+" => s"\\1.\\2"
julia> using CellListMap

julia> using PDBTools: read_pdb, coor

julia> positions = coor(read_pdb(CellListMap.argon_pdb_file));

julia> sys = ParticleSystem(
           xpositions = positions, 
           unitcell = [21.0, 21.0, 21.0],
           cutoff = 8.0, 
           output = 0.0, 
        );

julia> pairwise!((pair,output) -> output += pair.d2, sys)
43774.54367600001
```
## Two sets of particles

```jldoctest ;filter = r"(\\d*)\\.(\\d{2})\\d+" => s"\\1.\\2"
julia> using CellListMap, PDBTools

julia> xpositions = coor(read_pdb(CellListMap.argon_pdb_file))[1:50];

julia> ypositions = coor(read_pdb(CellListMap.argon_pdb_file))[51:100];

julia> sys = ParticleSystem(
           xpositions = xpositions, 
           ypositions = ypositions, 
           unitcell = [21.0, 21.0, 21.0],
           cutoff = 8.0, 
           output = 0.0, 
           parallel = false, # use true for parallelization
        );

julia> pairwise!((pair,output) -> output += pair.d2, sys)
21886.196785000004
```
"""
function ParticleSystem(;
        positions::SupportedCoordinatesTypes = nothing,
        xpositions::SupportedCoordinatesTypes = nothing,
        ypositions::SupportedCoordinatesTypes = nothing,
        unitcell::Union{AbstractVecOrMat, Nothing} = nothing,
        cutoff::Number,
        output::Any,
        output_name::Symbol = :output,
        parallel::Bool = true,
        nbatches::Tuple{Int, Int} = (0, 0),
        lcell = 1,
        validate_coordinates::Union{Nothing, Function} = _validate_coordinates,
    )
    # Set xpositions if positions was set
    if (isnothing(positions) && isnothing(xpositions)) || (!isnothing(positions) && !isnothing(xpositions))
        throw(ArgumentError("Either `positions` OR `xpositions` must be defined."))
    end
    xpositions = isnothing(positions) ? xpositions : positions
    # Check for simple input argument errors
    for input_array in (xpositions, ypositions)
        isnothing(input_array) && break
        if input_array isa AbstractMatrix
            dim = size(input_array, 1)
            if !(dim in (2, 3))
                throw(DimensionMismatch("Matrix of coordinates must have 2 or 3 rows, one for each dimension, got size: $(size(input_array))"))
            end
            input_array = reinterpret(reshape, SVector{dim, eltype(input_array)}, input_array)
        end
        if !isnothing(unitcell)
            DIM = if eltype(input_array) <: SVector
                length(eltype(input_array))
            else
                if length(input_array) == 0
                    # If the array is empty, we cannot determine the dimension, so we assume it is the same as the unit cell
                    size(unitcell, 1)
                else
                    length(input_array[1])
                end
            end
            if DIM != size(unitcell, 1)
                throw(DimensionMismatch("Dimension of the unit cell ($(size(unitcell, 1))) must match the dimension of the coordinates ($(length(eltype(input_array))))"))
            end
        end
    end
    # Single set of positions
    if isnothing(ypositions)
        unitcell = isnothing(unitcell) ? limits(xpositions; validate_coordinates) : unitcell
        _box = CellListMap.Box(unitcell, cutoff, lcell = lcell)
        _cell_list = CellListMap.CellList(xpositions, _box; parallel, nbatches, validate_coordinates)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list, :map)]
        output = _reset_all_output!(output, _output_threaded; reset = false)
        sys = ParticleSystem1{output_name}(xpositions, output, _box, _cell_list, _output_threaded, _aux, parallel, validate_coordinates)
        # Two sets of positions
    else
        unitcell = isnothing(unitcell) ? limits(xpositions, ypositions; validate_coordinates) : unitcell
        _box = CellListMap.Box(unitcell, cutoff, lcell = lcell)
        _cell_list = CellListMap.CellList(xpositions, ypositions, _box; parallel, nbatches, validate_coordinates)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list, :map)]
        output = _reset_all_output!(output, _output_threaded; reset = false)
        sys = ParticleSystem2{output_name}(xpositions, ypositions, output, _box, _cell_list, _output_threaded, _aux, parallel, validate_coordinates)
    end
    return sys
end

import Base: getproperty, propertynames
getproperty(sys::AbstractParticleSystem, s::Symbol) = getproperty(sys, Val(s))
getproperty(sys::AbstractParticleSystem, ::Val{S}) where {S} = getfield(sys, S)
# public properties
getproperty(sys::AbstractParticleSystem, ::Val{:unitcell}) = getfield(getfield(getfield(sys, :_box), :input_unit_cell), :matrix)
getproperty(sys::AbstractParticleSystem, ::Val{:cutoff}) = getfield(getfield(sys, :_box), :cutoff)
getproperty(sys::AbstractParticleSystem{OutputName}, ::Val{OutputName}) where {OutputName} = getfield(sys, :output)
propertynames(::AbstractParticleSystem{OutputName}) where {OutputName} =
    (:xpositions, :ypositions, :unitcell, :cutoff, :positions, :output, :parallel, OutputName)

import Base: setproperty!
# public properties
setproperty!(sys::AbstractParticleSystem, s::Symbol, x) = setproperty!(sys, Val(s), x)
setproperty!(sys::AbstractParticleSystem, ::Val{:unitcell}, x) = update_unitcell!(sys, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:cutoff}, x) = update_cutoff!(sys, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:parallel}, x) = setfield!(sys, :parallel, x)
# private properties
setproperty!(sys::AbstractParticleSystem, ::Val{:_box}, x) = setfield!(sys, :_box, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:_cell_list}, x) = setfield!(sys, :_cell_list, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:output}, x) = setfield!(sys, :output, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:_aux}, x) = setfield!(sys, :_aux, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:_output_threaded}, x) = setfield!(sys, :_output_threaded, x)

#=
    unitcelltype(sys::AbstractParticleSystem)

Returns the type of a unitcell from the `ParticleSystem` structure.

=#
CellListMap.unitcelltype(sys::AbstractParticleSystem) = unitcelltype(sys._box)

ParticleSystem1{OutputName}(v::V, o::O, b::B, c::C, vo::AbstractVector{O}, a::A, p::Bool, vc::VC) where {OutputName, V, O, B, C, A, VC} =
    ParticleSystem1{OutputName, V, O, B, C, A, VC}(v, o, b, c, vo, a, p, vc)
getproperty(sys::ParticleSystem1, ::Val{:positions}) = getfield(sys, :xpositions)
getproperty(sys::ParticleSystem2, ::Val{:positions}) = getfield(sys, :xpositions)

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem1{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "ParticleSystem1{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println(io)
    show(io_sub, mime, sys._cell_list)
    print(io, "\n    Parallelization auxiliary data set for $(nbatches(sys._cell_list, :build)) batch(es).")
    return print(io, "\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem2{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "ParticleSystem2{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println(io)
    show(io_sub, mime, sys._cell_list)
    print(io, "\n    Parallelization auxiliary data set for $(nbatches(sys._cell_list, :build)) batch(es).")
    return print(io, "\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

#
# Functions to copy, reset and reduce output variables, that must be implemented
# by the user for custom output types.
#
"""
    copy_output(x)

Defines how the `output` variable is copied. Identical to `Base.copy(x)`
and implemented for the types in `$(SupportedTypes)`.

Other custom output types must have their `copy_output` method implemented.

# Example

```julia
using CellListMap
# Custom data type
struct A x::Int end
# Custom output type (array of A)
output = [ A(0) for _ in 1:100 ]
# How to copy an array of `A`
CellListMap.copy_output(v::Vector{A}) = [ x for x in v ]

# Alternatively, in this case, one could have defined:
Base.copy(a::A) = a
CellListMap.copy_output(v::Vector{A}) = copy(v)
```

The user must guarantee that the copy is independent of the original array.
For many custom types it is possible to define 
```
CellListMap.copy_output(v::Vector{T}) where {T<:CustomType} = deepcopy(v)
```

"""
function copy_output(x)
    throw(
        ArgumentError(
            """\n
                No method matching `copy_output($(typeof(x)))`

                Please implement a method 
               
                CellListMap.copy_output(x::$(typeof(x)))

                with an appropriate way to copy the required output variable. Many times just
                defining `CellListMap.copy_output(x::$(typeof(x))) = deepcopy(x)` is ok. 
            """
        )
    )
end
copy_output(x::T) where {T <: SupportedTypes} = copy(x)
copy_output(x::AbstractVecOrMat{T}) where {T} = T[copy_output(el) for el in x]

"""
    reset_output(x)
    reset_output!(x)

Function that defines how to reset (or zero) the `output` variable. For `$(SupportedTypes)` it is 
implemented as `zero(x)`.

Other custom output types must have their `reset_output!` method implemented. 

The function *must* return the variable itself. If it is immutable,
a new instante of the variable must be created, with the reset value. 

!!! note
    By default, if
    `reset_output!` is defined for one element type, `reset_output!` is defined for arrays of that type
    by calling `reset_output!` for each element of the array.  The user must overload the `reset_output!` 
    function for the custom type array if that is not the desired behavior.

`reset_output` and `reset_output!` are aliases, and by convention `reset_output!` is preferred for mutable types.

# Example

In this example, we define a `reset_output` function that will set to `+Inf` the
minimum distance between particles (not always resetting means zeroing).

```jldoctest
julia> using CellListMap

julia> struct MinimumDistance d::Float64 end

julia> CellListMap.reset_output(x::MinimumDistance) = MinimumDistance(+Inf)

julia> x = MinimumDistance(1.0)
MinimumDistance(1.0)

julia> CellListMap.reset_output(x)
MinimumDistance(Inf)
```

See the `reducer` help entry for a complete example of how to use `reset_output`.

"""
function reset_output!(x)
    throw(
        ArgumentError(
            """\n
                No method matching `reset_output!($(typeof(x)))`

                Please add a method 
                
                CellListMap.reset_output!(x::$(typeof(x)))
                
                with the appropriate way to reset (zero) the data of the output variables.

                The reset_output! methods **must** return the output variable to
                conform with the interface, even if the variable is mutable. 
            """
        )
    )
end
reset_output!(x::T) where {T <: SupportedTypes} = zero(x)
function reset_output!(x::AbstractVecOrMat{T}) where {T}
    for i in eachindex(x)
        x[i] = reset_output!(x[i])
    end
    return x
end
const reset_output = reset_output!
#=
    _reset_all_output!(output, output_threaded)

Function that resets the output variable and the threaded copies of it.

=#
function _reset_all_output!(output, output_threaded; reset::Bool)
    if reset
        output = reset_output!(output)
    end
    for i in eachindex(output_threaded)
        output_threaded[i] = reset_output!(output_threaded[i])
    end
    return output
end
"""
    reducer(x,y)
    reducer!(x,y)

Defines how to reduce (combine, or merge) to variables computed in parallel
to obtain a single instance of the variable with the reduced result. 

`reducer` and `reducer!` are aliases, and `reducer!` is preferred, by convention
for mutating functions.

The most common `reducer` is the sum, and this is how it is implemented for
`$(SupportedTypes)`. For example, when computing energies, or forces,
the total energy is the sum of the energies. The force on one particle is the sum of the
forces between the particle and every other particle. Thus, the implemented reducer is
the sum: 

```
reducer(x,y) = +(x,y)
```

However, in  many cases, reduction must be done differently. For instance, if the minimum
distance between particles is to be computed, it is interesting to define a custom type
and associated reducer. For example:

```
struct MinimumDistance d::Float64 end
reducer(x::MinimumDistance, y::MinimumDistance) = MinimumDistance(min(x.d, y.d))
```

The overloading of `reducer` allows the use of parallel computations for custom, 
complex data types, containing different types of variables, fields, or sizes.

The appropriate behavior of the reducer should be carefully inspected by the user
to avoid spurious results. 

# Example

In this example we show how to obtain the minimum distance among argon atoms
in a simulation box.

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
julia> using CellListMap, PDBTools

julia> positions = coor(read_pdb(CellListMap.argon_pdb_file));

julia> struct MinimumDistance d::Float64 end # Custom output type

julia> CellListMap.copy_output(d::MinimumDistance) = MinimumDistance(d.d) # Custom copy function for `Out`

julia> CellListMap.reset_output(d::MinimumDistance) = MinimumDistance(+Inf) # How to reset an array with elements of type `MinimumDistance`

julia> CellListMap.reducer(md1::MinimumDistance, md2::MinimumDistance) = MinimumDistance(min(md1.d, md2.d)) # Custom reduction function

julia> # Construct the system
       sys = ParticleSystem(;
           positions = positions,
           unitcell = [21,21,21],
           cutoff = 8.0,
           output = MinimumDistance(+Inf),
       );

julia> # Obtain the minimum distance between atoms:
       pairwise!((pair,output) -> pair.d < output.d ? MinimumDistance(pair.d) : output, sys)
MinimumDistance(2.1991993997816563)
```

"""
function reducer!(x, y)
    throw(
        ArgumentError(
            """\n
                No method matching `reducer!($(typeof(x)),$(typeof(y)))`

                Please implement a method 
                
                CellListMap.reducer(x::$(typeof(x)),y::$(typeof(y)))
                
                with the appropriate way to combine two instances of the type (summing, keeping
                the minimum, etc), such that threaded computations can be reduced.

            """
        )
    )
end
reducer!(x::T, y::T) where {T <: SupportedTypes} = +(x, y)
const reducer = reducer!

#=
    reduce_output!(reducer::Function, output, output_threaded)

# Extended help

Function that defines how to reduce the vector of `output` variables, after a threaded
computation. This function is implemented for `output` variables that are numbers, 
and vectors or arrays of number of static arrays, as the sum of the values of the 
threaded computations, which is the most common application, found in computing
forces, energies, etc. 

It may be interesting to implement custom `CellListMap.reduce_output!` function for other types 
of output variables, considering:

- The arguments of the function must be the return `output` value and a vector 
  `output_threaded` of `output` variables, which is created (automatically) by copying
  the output the number of times necessary for the multi-threaded computation. 

- The function *must* return the `output` variable, independently of it being mutable
  or immutable.

`reduce_output` is an alias to `reduce_output!` that can be used for consistency if the `output`
variable is immutable.

```

=#
function reduce_output!(reducer::Function, output::T, output_threaded::Vector{T}) where {T}
    #output = reset_output!(output)
    for ibatch in eachindex(output_threaded)
        output = reducer(output, output_threaded[ibatch])
    end
    return output
end
function reduce_output!(
        reducer::Function,
        output::AbstractVecOrMat{T},
        output_threaded::Vector{<:AbstractVecOrMat{T}}
    ) where {T}
    output = reset_output!(output)
    for ibatch in eachindex(output_threaded)
        for i in eachindex(output, output_threaded[ibatch])
            output[i] = reducer(output[i], output_threaded[ibatch][i])
        end
    end
    return output
end

"""
    resize_output!(sys::AbstractParticleSystem, n::Int)

Resizes the output array and the auxiliary output arrays used
for multithreading, if the number of particles of the system changed.

This function must be implemented by the user if the output variable is a 
vector whose length is dependent on the number of particles. For example,
if the output is a vector of forces acting on each particle, the output
vector must be resized if the number of particles changes. 

This function *must* be used in that case, to guarantee that the 
auxiliary arrays used for multi-threading are resized accordingly. 

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

#
# This interface is needed to generate random particle coordinates in ComplexMixtures.jl
#
#=
    get_computing_box(sys::AbstractParticleSystem)

Retrieves the computing box of the system. The computing box is large enough to
contain all coordinates of the particles, plus the cutoff.

=#
get_computing_box(sys::AbstractParticleSystem) = sys._box.computing_box

#=
    UpdateParticleSystem!

Updates the cell lists for periodic systems.

=#
function UpdateParticleSystem!(sys::ParticleSystem1, update_lists::Bool = true)
    if update_lists
        if unitcelltype(sys) == NonPeriodicCell
            sys._box = Box(limits(sys.xpositions), sys.cutoff)
        end
        n_particles_changed = length(sys.xpositions) != sys._cell_list.n_real_particles
        sys._cell_list = CellListMap.UpdateCellList!(
            sys.xpositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel = sys.parallel,
            validate_coordinates = sys.validate_coordinates,
        )
        if n_particles_changed
            _old_nbatches = nbatches(sys)
            sys._cell_list = update_number_of_batches!(sys._cell_list; parallel = sys.parallel)
            _new_nbatches = nbatches(sys)
            if _old_nbatches != _new_nbatches
                sys._aux = CellListMap.AuxThreaded(sys._cell_list)
                sys._output_threaded = [copy_output(sys.output) for _ in 1:_new_nbatches[2]]
            end
        end
    end
    return sys
end

function UpdateParticleSystem!(sys::ParticleSystem2, update_lists::Bool = true)
    if update_lists
        if unitcelltype(sys) == NonPeriodicCell
            sys._box = Box(limits(sys.xpositions, sys.ypositions), sys.cutoff)
        end
        n_particles_changed = (min(length(sys.xpositions), length(sys.ypositions)) != sys._cell_list.small_set.n_real_particles) ||
            (max(length(sys.xpositions), length(sys.ypositions)) != sys._cell_list.large_set.n_real_particles)
        sys._cell_list = CellListMap.UpdateCellList!(
            sys.xpositions,
            sys.ypositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel = sys.parallel,
            validate_coordinates = sys.validate_coordinates,
        )
        if n_particles_changed
            _old_nbatches = nbatches(sys)
            sys._cell_list = update_number_of_batches!(sys._cell_list; parallel = sys.parallel)
            _new_nbatches = nbatches(sys)
            if _old_nbatches != _new_nbatches
                sys._aux = CellListMap.AuxThreaded(sys._cell_list)
                sys._output_threaded = [copy_output(sys.output) for _ in 1:_new_nbatches[2]]
            end
        end
    end
    return sys
end

# Return the number of batches for ParticleSystems
nbatches(sys::ParticleSystem1) = nbatches(sys._cell_list)
nbatches(sys::ParticleSystem2) = nbatches(sys._cell_list.small_set)

"""
    pairwise!(
        f::Function, system::AbstractParticleSystem; 
        show_progress=true, update_lists=true, reset=true,
    )

Function that maps the `f` function into all pairs of particles of
`system` that are found to be within the `cutoff`. 

The function `f` receives a `NeighborPair` struct and the output:
```
function f(pair, output)
    # pair.i, pair.j: indices of the particles
    # pair.x, pair.y: coordinates (minimum-image adjusted)
    # pair.d: distance between particles
    # pair.d2: squared distance
    # update output
    return output
end
```

Thread-safety is taken care automatically in parallel executions.

`pairwise` is an alias to `pairwise!` for syntax consistency
when the `output` variable is immutable.

If `update_lists` is `false`, the cell lists will not be recomputed,
this may be useful for computing a different function from the same
coordinates.

If `reset` is set to `false`, the value of `system.output` will not be
set to `zero(typeof(system.output))` before the new accumulation.

# Example

In this example we compute the sum of `1/(1+d)` where `d` is the
distance between particles of a set, for `d < cutoff`. 

```julia-repl
julia> sys = ParticleSystem(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           output = 0.0
           );

julia> pairwise!((pair, output) -> output += 1 / (1 + pair.d), sys)
1870.0274887950268
```

"""
function pairwise!(
        f::F,
        sys::AbstractParticleSystem;
        update_lists::Bool = true,
        show_progress::Bool = false,
        reset::Bool = true,
    ) where {F <: Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded; reset)
    UpdateParticleSystem!(sys, update_lists)
    sys.output = CellListMap.pairwise!(
        f, sys.output, sys._box, sys._cell_list;
        output_threaded = sys._output_threaded,
        parallel = sys.parallel,
        reduce = (output, output_threaded) -> reduce_output!(reducer, output, output_threaded),
        show_progress = show_progress
    )
    return sys.output
end

#
# Cross-computations when only one cell list was computed
#
"""
    pairwise!(f::Function, x::AbstractVector{<:AbstractVector}, sys::ParticleSystem1; kargs...)
    pairwise!(f::Function, x::AbstractMatrix, sys::ParticleSystem1; kargs...)

Evaluate function f for pairs in two independent sets of particles, where the `sys::ParticleSystem1` object
contains the previously computed cell lists of one set of particles, and the second set is given by the
array of positions `x`.

This function can be advantageous over computing the interactions with `CellListPair`, because here the
cell lists are only computed for one set. This is advantageous in two situations:

    1. The second set of particles is not changing, and the first set is changing. Thus, the cell lists
       of the second set can be computed only once.
    2. One of the sets is much smaller than the other. In this case, computing the cell lists of the largest
       set might be too expensive. Construct the `ParticleSystem` object for the smallest set, and use this
       function to compute the interactions with the largest set.

## Keyword arguments:

- `show_progress::Bool=false`: Show progress bar.
- `update_lists::Bool=true`: Update the cell lists or not. If the positions of the `ParticleSystem1` object
   have not changed, it is not necessary to update the cell lists.
- `reset::Bool=true`: If set to `false` the value of `sys.output` will not be set to `zero(typeof(sys.output)`,
   and the result will be accumulated

## Example

```julia-repl
julia> using CellListMap, StaticArrays

julia> x = rand(SVector{3,Float64}, 1000);

julia> sys = ParticleSystem(positions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0);

julia> y = rand(SVector{3,Float64}, 100);

julia> pairwise!((pair, output) -> output + pair.d, y, sys; update_lists=false) # Compute the sum of the distances of x and y
31.121496300032163

julia> z = rand(SVector{3,Float64}, 200);

julia> pairwise!((pair, output) -> output + pair.d, z, sys; update_lists=false) # Compute the sum of the distances x and z
63.57860511891242
```

"""
function pairwise!(
        f::F, x::AbstractVecOrMat, sys::ParticleSystem1;
        show_progress::Bool = false, update_lists::Bool = true, reset::Bool = true,
    ) where {F <: Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded; reset)
    UpdateParticleSystem!(sys, update_lists)
    sys.output = pairwise!(
        f, sys.output, sys._box, x, sys._cell_list;
        output_threaded = sys._output_threaded,
        reduce = (output, output_threaded) -> reduce_output!(reducer, output, output_threaded),
        parallel = sys.parallel, show_progress,
    )
    return sys.output
end
