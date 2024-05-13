export ParticleSystem
export update_cutoff!
export update_unitcell!
export resize_output!
@compat public copy_output, reset_output!, reset_output, reducer, reducer!

# Types of variables that have support for multi-threading without having 
# to explicit add methods to copy_output, reset_output!, and reducer functions.
const SupportedTypes = Union{Number,SVector,FieldVector}

# Supported types for coordinates
const SupportedCoordinatesTypes = Union{Nothing,AbstractVector{<:AbstractVector},AbstractMatrix}

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
        autoswap::Bool = true
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
full unit cell matrix for `Triclinic` cells), the cutoff used for the
construction of the cell lists and the output variable of the calculations.
If unitcell == nothing, the system is considered not-periodic, in which
case artificial periodic boundaries will be built such that images 
are farther from each other than the cutoff.

`output_name` can be set to a symbol that best identifies the output variable.
For instance, if `output_name=:forces`, the forces can be retrieved from the
structure using the `system.forces` notation.

The `parallel` and `nbatches` flags control the parallelization scheme of
computations (see https://m3g.github.io/CellListMap.jl/stable/parallelization/#Number-of-batches)).
By default the parallelization is turned on and `nbatches` is set with heuristics
that may provide good efficiency in most cases. `autoswap = false` will guarantee that
the cell lists will be buitl for the `ypositions` (by default they are constructed
for the smallest set, which is faster).

# Example

In these examples, we compute the sum of the squared distances between
the particles that are within the cutoff:

## Single set of particles

```jldoctest ;filter = r"(\\d*)\\.(\\d)\\d+" => s"\\1.\\2"
julia> using CellListMap

julia> using PDBTools: readPDB, coor

julia> positions = coor(readPDB(CellListMap.argon_pdb_file));

julia> sys = ParticleSystem(
           xpositions = positions, 
           unitcell = [21.0, 21.0, 21.0],
           cutoff = 8.0, 
           output = 0.0, 
        );

julia> map_pairwise!((x,y,i,j,d2,output) -> output += d2, sys)
43774.54367600001
```
## Two sets of particles

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
julia> using CellListMap, PDBTools

julia> xpositions = coor(readPDB(CellListMap.argon_pdb_file))[1:50];

julia> ypositions = coor(readPDB(CellListMap.argon_pdb_file))[51:100];

julia> sys = ParticleSystem(
           xpositions = xpositions, 
           ypositions = ypositions, 
           unitcell = [21.0, 21.0, 21.0],
           cutoff = 8.0, 
           output = 0.0, 
           parallel = false, # use true for parallelization
        );

julia> map_pairwise!((x,y,i,j,d2,output) -> output += d2, sys)
21886.196785000004
```
"""
function ParticleSystem(;
    positions::SupportedCoordinatesTypes=nothing,
    xpositions::SupportedCoordinatesTypes=nothing,
    ypositions::SupportedCoordinatesTypes=nothing,
    unitcell::Union{AbstractVecOrMat,Nothing}=nothing,
    cutoff::Number,
    output::Any,
    output_name::Symbol=:output,
    parallel::Bool=true,
    nbatches::Tuple{Int,Int}=(0, 0),
    lcell=1,
    autoswap::Bool=true
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
            input_array = reinterpret(reshape, SVector{dim,eltype(input_array)}, input_array)
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
        unitcell = isnothing(unitcell) ? limits(xpositions) : unitcell
        _box = CellListMap.Box(unitcell, cutoff, lcell=lcell)
        _cell_list = CellListMap.CellList(xpositions, _box; parallel=parallel, nbatches=nbatches)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = ParticleSystem1{output_name}(xpositions, output, _box, _cell_list, _output_threaded, _aux, parallel)
        # Two sets of positions
    else
        unitcell = isnothing(unitcell) ? limits(xpositions, ypositions) : unitcell
        _box = CellListMap.Box(unitcell, cutoff, lcell=lcell)
        _cell_list = CellListMap.CellList(xpositions, ypositions, _box; parallel=parallel, nbatches=nbatches, autoswap=autoswap)
        _aux = CellListMap.AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list)]
        output = _reset_all_output!(output, _output_threaded)
        sys = ParticleSystem2{output_name}(xpositions, ypositions, output, _box, _cell_list, _output_threaded, _aux, parallel)
    end
    return sys
end

# Abstract type only for cleaner dispatch
abstract type AbstractParticleSystem{OutputName} end

import Base: getproperty, propertynames
getproperty(sys::AbstractParticleSystem, s::Symbol) = getproperty(sys, Val(s))
getproperty(sys::AbstractParticleSystem, s::Val{S}) where {S} = getfield(sys, S)
# public properties
getproperty(sys::AbstractParticleSystem, ::Val{:unitcell}) = getfield(getfield(getfield(sys, :_box), :input_unit_cell), :matrix)
getproperty(sys::AbstractParticleSystem, ::Val{:cutoff}) = getfield(getfield(sys, :_box), :cutoff)
getproperty(sys::AbstractParticleSystem{OutputName}, ::Val{OutputName}) where {OutputName} = getfield(sys, :output)
propertynames(sys::AbstractParticleSystem{OutputName}) where {OutputName} =
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

"""
    unitcelltype(sys::AbstractParticleSystem)

Returns the type of a unitcell from the `ParticleSystem` structure.

"""
CellListMap.unitcelltype(sys::AbstractParticleSystem) = unitcelltype(sys._box)

@testitem "ParticleSystem properties" begin

    using CellListMap
    using StaticArrays
    sys = ParticleSystem(
        positions=rand(SVector{3,Float64}, 1000),
        cutoff=0.1,
        unitcell=[1, 1, 1],
        output=0.0,
        output_name=:test
    )
    @test length(sys.positions) == 1000
    @test sys.cutoff == 0.1
    @test sys.unitcell == @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test sys.output == 0
    @test sys.test == 0
    @test sys.parallel == true

    sys.parallel = false
    @test sys.parallel == false
    sys.cutoff = 0.2
    @test sys.cutoff == 0.2
    sys.positions[1] = SVector(0.0, 0.0, 0.0)
    @test sys.positions[1] == SVector(0.0, 0.0, 0.0)
    sys.unitcell = [1.2, 1.2, 1.2]
    @test sys.unitcell == @SMatrix [1.2 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.2]

    # test the construction with pathologically few particles
    for x in [
        SVector{3,Float64}[],
        Vector{Float64}[],
        Matrix{Float64}(undef, 3, 0),
        [rand(SVector{3,Float64})],
        [rand(3)],
        rand(3, 1)
    ]
        _sys = ParticleSystem(
            positions=x,
            cutoff=0.1,
            unitcell=[1, 1, 1],
            output=0.0,
            output_name=:test
        )
        @test CellListMap.map_pairwise((x, y, i, j, d2, out) -> out += d2, _sys) == 0.0
    end

    # unitcell type
    x = rand(SVector{3,Float64}, 100)
    @test unitcelltype(ParticleSystem(positions=x, cutoff=0.1, unitcell=[1, 1, 1], output=0.0)) == OrthorhombicCell
    @test unitcelltype(ParticleSystem(positions=x, cutoff=0.1, unitcell=[1 0 0; 0 1 0; 0 0 1], output=0.0)) == TriclinicCell
    @test unitcelltype(ParticleSystem(positions=x, cutoff=0.1, output=0.0)) == NonPeriodicCell

    # Argument errors
    @test_throws ArgumentError ParticleSystem(
        positions=rand(SVector{3,Float64}, 100),
        xpositions=rand(SVector{3,Float64}, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws ArgumentError ParticleSystem(
        ypositions=rand(SVector{3,Float64}, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        positions=rand(1, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        xpositions=rand(1, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        xpositions=rand(2, 100),
        ypositions=rand(1, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        positions=rand(2, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        xpositions=rand(2, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        xpositions=rand(2, 100),
        ypositions=rand(2, 100),
        cutoff=0.1, unitcell=[1, 1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        positions=rand(3, 100),
        cutoff=0.1, unitcell=[1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        xpositions=rand(3, 100),
        cutoff=0.1, unitcell=[1, 1], output=0.0,
    )
    @test_throws DimensionMismatch ParticleSystem(
        xpositions=rand(3, 100),
        ypositions=rand(3, 100),
        cutoff=0.1, unitcell=[1, 1], output=0.0,
    )

end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `map_pairwise!` computations,
for systems with one set of positions (thus, replacing the loops over `N(N-1)` 
pairs of particles of the set). 

The `xpositions`, `output`, and `parallel` fields are considered part of the API,
and you can retrive or mutate `xpositions`, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fileds of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `ParticleSystem1` structure
is done through the `ParticleSystem(;xpositions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct ParticleSystem1{OutputName,V,O,B,C,A} <: AbstractParticleSystem{OutputName}
    xpositions::V
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
end
ParticleSystem1{OutputName}(v::V, o::O, b::B, c::C, vo::AbstractVector{O}, a::A, p::Bool) where {OutputName,V,O,B,C,A} =
    ParticleSystem1{OutputName,V,O,B,C,A}(v, o, b, c, vo, a, p)
getproperty(sys::ParticleSystem1, ::Val{:positions}) = getfield(sys, :xpositions)

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Structure that carries the information necessary for `map_pairwise!` computations,
for systems with two set of positions (thus, replacing the loops over `N×M` 
pairs of particles, being `N` and `M` the number of particles of each set).

The `xpositions`, `ypositions`, `output`, and `parallel` fields are considered part of the API,
and you can retrive or mutate positions, retrieve the `output` or its elements,
and set the computation to use or not parallelization by directly accessing these
elements.

The other fileds of the structure (starting with `_`) are internal and must not 
be modified or accessed directly. The construction of the `ParticleSystem1` structure
is done through the `ParticleSystem(;xpositions, ypositions, unitcell, cutoff, output)` 
auxiliary function.

"""
mutable struct ParticleSystem2{OutputName,V,O,B,C,A} <: AbstractParticleSystem{OutputName}
    xpositions::V
    ypositions::V
    output::O
    _box::B
    _cell_list::C
    _output_threaded::Vector{O}
    _aux::A
    parallel::Bool
end
ParticleSystem2{OutputName}(vx::V, vy::V, o::O, b::B, c::C, vo::Vector{O}, a::A, p::Bool) where {OutputName,V,O,B,C,A} =
    ParticleSystem2{OutputName,V,O,B,C,A}(vx, vy, o, b, c, vo, a, p)

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem1{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "ParticleSystem1{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println(io)
    show(io_sub, mime, sys._cell_list)
    println(io, "\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.nbatches)
    print(io, "\n    Type of output variable ($OutputName): $(typeof(sys.output))")
end

function Base.show(io::IO, mime::MIME"text/plain", sys::ParticleSystem2{OutputName}) where {OutputName}
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 4)
    N = size(sys.unitcell, 1)
    println(io, "ParticleSystem2{$OutputName} of dimension $N, composed of:")
    show(IOContext(io, :indent => indent + 4), mime, sys._box)
    println(io)
    show(io_sub, mime, sys._cell_list)
    println(io, "\n    Parallelization auxiliary data set for: ")
    show(io_sub, mime, sys._cell_list.target.nbatches)
    print(io, "\n    Type of output variable ($OutputName): $(typeof(sys.output))")
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

# Alternativelly, in this case, one could have defined:
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
    throw(ArgumentError("""\n
        No method matching `copy_output($(typeof(x)))`

        Please implement a method 
       
        CellListMap.copy_output(x::$(typeof(x)))

        with an appropriate way to copy the required output variable. Many times just
        defining `output_copy(x::$(typeof(x))) = deepcopy(x)` is ok. 
    """))
end
copy_output(x::T) where {T<:SupportedTypes} = copy(x)
copy_output(x::AbstractVecOrMat{T}) where {T} = T[copy_output(el) for el in x]

"""
    reset_output(x)
    reset_output!(x)

Function that defines how to reset (or zero) the `output` variable. For `$(SupportedTypes)` it is 
implemented as `zero(x)`, and for `AbstractVecOrMat` containers of `Number`s or `SVector`s
it is implemented as `fill!(x, zero(eltype(x))`.

Other custom output types must have their `reset_output!` method implemented.

If the variable is mutable, the function *must* return the variable itself. If it is immutable,
a new instante of the variable must be created, with the reset value. 

`reset_output` and `reset_output!` are aliases, and `reset_output!` is preferred, by convention
for mutating functions.

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
    throw(ArgumentError("""\n
        No method matching `reset_output!($(typeof(x)))`

        Please add a method 
        
        CellListMap.reset_output!(x::$(typeof(x)))
        
        with the appropriate way to reset (zero) the data of the output variables.

        The reset_output! methods **must** return the output variable to
        conform with the interface, even if the variable is mutable. 
    """))
end
reset_output!(x::T) where {T<:SupportedTypes} = zero(x)
reset_output!(x::AbstractVecOrMat{T}) where {T} = fill!(x, reset_output!(x[begin]))
const reset_output = reset_output!

#=
    _reset_all_output!(output, output_threaded)

Function that resets the output variable and the threaded copies of it.

=#
function _reset_all_output!(output, output_threaded)
    output = reset_output!(output)
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

The most commont `reducer` is the sum, and this is how it is implemented for
`$(SupportedTypes)`. For example, when computin energies, or forces,
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

The appropriate behavior of the reducer should be carefuly inspected by the user
to avoid spurious results. 

# Example

In this example we show how to obtain the minimum distance among argon atoms
in a simulation box.

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
julia> using CellListMap, PDBTools

julia> positions = coor(readPDB(CellListMap.argon_pdb_file));

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
       map_pairwise!((x,y,i,j,d2,output) -> sqrt(d2) < output.d ? MinimumDistance(sqrt(d2)) : output, sys)
MinimumDistance(2.1991993997816563)
```

"""
function reducer!(x, y)
    throw(ArgumentError("""\n
        No method matching `reducer!($(typeof(x)),$(typeof(y)))`

        Please implement a method 
        
        CellListMap.reducer(x::$(typeof(x)),y::$(typeof(y)))
        
        with the appropriate way to combine two instances of the type (summing, keeping
        the minimum, etc), such that threaded computations can be reduced.
    """))
end
reducer!(x::T, y::T) where {T<:SupportedTypes} = +(x, y)
const reducer = reducer!

@testitem "reducer method basics and errors" begin
    using CellListMap
    @test CellListMap.reducer(1, 2) == 3
    @test CellListMap.copy_output(1) == 1
    @test CellListMap.reset_output(1) == 0
    struct A end
    @test_throws ArgumentError CellListMap.reset_output(A())
    @test_throws ArgumentError CellListMap.copy_output(A())
    @test_throws ArgumentError CellListMap.reducer(A(), A())
end

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
    output = reset_output!(output)
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
`N` is the dimension of the sytem (2D or 3D).

This function can be used to update the system geometry in iterative schemes,
where the size of the simulation box changes during the simulation.

!!! note
    Manual updating of the unit cell of non-periodic systems is not allowed.

# Example

```jldoctest ;filter = r"batches.*" => ""  
julia> using CellListMap, StaticArrays, PDBTools

julia> xpositions = coor(readPDB(CellListMap.argon_pdb_file));

julia> sys = ParticleSystem(
           xpositions = xpositions,
           unitcell=[21,21,21], 
           cutoff = 8.0, 
           output = 0.0
       );

julia> update_unitcell!(sys, [30.0, 30.0, 30.0])
ParticleSystem1{output} of dimension 3, composed of:
    Box{OrthorhombicCell, 3}
      unit cell matrix = [ 30.0 0.0 0.0; 0.0 30.0 0.0; 0.0 0.0 30.0 ]
      cutoff = 8.0
      number of computing cells on each dimension = [6, 6, 6]
      computing cell sizes = [10.0, 10.0, 10.0] (lcell: 1)
      Total number of cells = 216
    CellList{3, Float64}
      100 real particles.
      8 cells with real particles.
      800 particles in computing box, including images.
    Parallelization auxiliary data set for:
      Number of batches for cell list construction: 1
      Number of batches for function mapping: 1
    Type of output variable (output): Float64

```

"""
function update_unitcell!(sys, unitcell)
    if unitcelltype(sys) == NonPeriodicCell
        throw(ArgumentError("""\n
            Manual updating of the unit cell of non-periodic systems is not allowed.
        """))
    end
    sys._box = update_box(sys._box; unitcell=unitcell)
    return sys
end

@testitem "update_unitcell!" begin
    using BenchmarkTools
    using LinearAlgebra: diag
    using StaticArrays
    using CellListMap
    x = rand(SVector{3,Float64}, 1000)
    sys1 = ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_unitcell!(sys1, SVector(2, 2, 2))
    @test diag(sys1.unitcell) == [2, 2, 2]
    a = @ballocated update_unitcell!($sys1, SVector(2, 2, 2)) evals = 1 samples = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys2 = ParticleSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_unitcell!(sys2, SVector(2, 2, 2))
    @test diag(sys2.unitcell) == [2, 2, 2]
    a = @ballocated update_unitcell!($sys2, SVector(2, 2, 2)) evals = 1 samples = 1
    @test a == 0
    # Test throwing error on updating non-periodic unit cells
    sys = ParticleSystem(xpositions=x, cutoff=0.1, output=0.0)
    @test_throws ArgumentError update_unitcell!(sys, [1, 1, 1])
    sys = ParticleSystem(xpositions=x, ypositions=y, cutoff=0.1, output=0.0)
    @test_throws ArgumentError update_unitcell!(sys, [1, 1, 1])
end

"""
    update_cutoff!(system, cutoff)

Function to update the `cutoff`` of the system. 

This function can be used to update the system geometry in iterative schemes.

# Example

Here we initialize a particle system with a cutoff of `8.0` and then update
the cutoff to `10.0`. 

```jldoctest ; filter = r"batches.*" => ""
julia> using CellListMap, PDBTools

julia> x = coor(readPDB(CellListMap.argon_pdb_file));

julia> sys = ParticleSystem(
           xpositions = x, 
           unitcell=[21.0,21.0,21.0], 
           cutoff = 8.0, 
           output = 0.0
       );

julia> update_cutoff!(sys, 10.0)
ParticleSystem1{output} of dimension 3, composed of:
    Box{OrthorhombicCell, 3}
      unit cell matrix = [ 21.0 0.0 0.0; 0.0 21.0 0.0; 0.0 0.0 21.0 ]
      cutoff = 10.0
      number of computing cells on each dimension = [5, 5, 5]
      computing cell sizes = [10.5, 10.5, 10.5] (lcell: 1)
      Total number of cells = 125
    CellList{3, Float64}
      100 real particles.
      8 cells with real particles.
      800 particles in computing box, including images.
    Parallelization auxiliary data set for:
      Number of batches for cell list construction: 8
      Number of batches for function mapping: 8
    Type of output variable (output): Float64
```
"""
function update_cutoff!(sys::ParticleSystem1, cutoff)
    if unitcelltype(sys) == NonPeriodicCell
        sys._box = Box(limits(sys.xpositions), cutoff)
    end
    sys._box = update_box(sys._box; cutoff=cutoff)
    return sys
end
function update_cutoff!(sys::ParticleSystem2, cutoff)
    if unitcelltype(sys) == NonPeriodicCell
        sys._box = Box(limits(sys.xpositions, sys.ypositions), cutoff)
    end
    sys._box = update_box(sys._box; cutoff=cutoff)
    return sys
end

@testitem "update_cutoff!" begin
    using BenchmarkTools
    using StaticArrays
    using CellListMap
    using PDBTools
    x = rand(SVector{3,Float64}, 1000)
    sys1 = ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_cutoff!(sys1, 0.2)
    @test sys1.cutoff == 0.2
    a = @ballocated update_cutoff!($sys1, 0.1) evals = 1 samples = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys2 = ParticleSystem(xpositions=x, ypositions=y, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    update_cutoff!(sys2, 0.2)
    @test sys2.cutoff == 0.2
    a = @ballocated update_cutoff!($sys2, 0.1) evals = 1 samples = 1
    @test a == 0

    # Update cutoff of non-periodic systems
    x = coor(readPDB(CellListMap.argon_pdb_file))
    sys1 = ParticleSystem(xpositions=x, cutoff=8.0, output=0.0)
    @test unitcelltype(sys1) == NonPeriodicCell
    @test sys1.unitcell ≈ [ 35.63 0.0 0.0; 0.0 35.76 0.0; 0.0 0.0 35.79 ] atol = 1e-2
    update_cutoff!(sys1, 10.0)
    @test sys1.unitcell ≈ [ 39.83 0.0 0.0; 0.0 39.96 0.0; 0.0 0.0 39.99 ] atol = 1e-2
    a = @ballocated update_cutoff!($sys1, 8.0) evals = 1 samples = 1
    @test a == 0
    sys2 = ParticleSystem(xpositions=x[1:50], ypositions=x[51:100], cutoff=8.0, output=0.0)
    @test unitcelltype(sys2) == NonPeriodicCell
    @test sys2.unitcell ≈ [ 35.63 0.0 0.0; 0.0 35.76 0.0; 0.0 0.0 35.79 ] atol = 1e-2
    update_cutoff!(sys2, 10.0)
    @test sys2.unitcell ≈ [ 39.83 0.0 0.0; 0.0 39.96 0.0; 0.0 0.0 39.99 ] atol = 1e-2
    a = @ballocated update_cutoff!($sys2, 8.0) evals = 1 samples = 1
    @test a == 0
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
@testitem "get_computing_box" begin
    using StaticArrays
    using CellListMap
    x = rand(SVector{3,Float64}, 1000)
    sys = ParticleSystem(xpositions=x, unitcell=[1, 1, 1], cutoff=0.1, output=0.0)
    @test CellListMap.get_computing_box(sys) == ([-0.1, -0.1, -0.1], [1.1, 1.1, 1.1])
end

#=
    UpdateParticleSystem!

Updates the cell lists for periodic systems.

=#
function UpdateParticleSystem!(sys::ParticleSystem1, update_lists::Bool=true)
    if update_lists
        if unitcelltype(sys) == NonPeriodicCell
            sys._box = Box(limits(sys.xpositions), sys.cutoff)
        end
        sys._cell_list = CellListMap.UpdateCellList!(
            sys.xpositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel=sys.parallel
        )
    end
    return sys
end

function _update_ref_positions!(cl::CellListPair{V,N,T,Swap}, sys) where {V,N,T,Swap<:NotSwapped}
    resize!(cl.ref, length(sys.xpositions))
    cl.ref .= sys.xpositions
end
function _update_ref_positions!(::CellListPair{V,N,T,Swap}, sys) where {V,N,T,Swap<:Swapped}
    throw(ArgumentError("update_lists == false requires autoswap == false for 2-set systems."))
end
function UpdateParticleSystem!(sys::ParticleSystem2, update_lists::Bool=true)
    if update_lists
        if unitcelltype(sys) == NonPeriodicCell
            sys._box = Box(limits(sys.xpositions, sys.ypositions), sys.cutoff)
        end
        sys._cell_list = CellListMap.UpdateCellList!(
            sys.xpositions,
            sys.ypositions,
            sys._box,
            sys._cell_list,
            sys._aux;
            parallel=sys.parallel
        )
    else
        # Always update the reference set positions (the cell lists of the target set are not updated)
        _update_ref_positions!(sys._cell_list, sys)
    end
    return sys
end

# this updates must be non-allocating in the serial case
@testitem "UpdateParticleSystem!" begin
    using BenchmarkTools
    using StaticArrays
    using CellListMap
    x = rand(SVector{3,Float64}, 1000)
    sys = ParticleSystem(xpositions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys = ParticleSystem(xpositions=x, ypositions=y, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0

    # Test construction with more general abstract vectors
    x = @view(x[1:500])
    sys = ParticleSystem(xpositions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0
    y = @view(y[1:500])
    sys = ParticleSystem(xpositions=x, ypositions=y, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0

    # Update with matrices
    x = rand(3, 500)
    sys = ParticleSystem(xpositions=x, unitcell=[1.0, 1.0, 1.0], cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0

    # Update non-periodic system
    x = rand(SVector{3,Float64}, 1000)
    sys = ParticleSystem(xpositions=x, cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0
    y = rand(SVector{3,Float64}, 1000)
    sys = ParticleSystem(xpositions=x, ypositions=y, cutoff=0.1, output=0.0, parallel=false)
    a = @ballocated CellListMap.UpdateParticleSystem!($sys) samples = 1 evals = 1
    @test a == 0

    # Throw error when trying to *not* update lists with autoswap on:
    sys = ParticleSystem(
        xpositions = rand(SVector{2,Float64}, 100),
        ypositions = rand(SVector{2,Float64}, 200),
        unitcell = [1,1],
        cutoff = 0.1,
        output=0.0,
        autoswap=true,
    )
    @test_throws ArgumentError map_pairwise!((_, _, _, _, d2, u) -> u += d2, sys, update_lists=false)

end

"""
    map_pairwise!(
        f::Function, system::AbstractParticleSystem; 
        show_progress = true, update_lists = true
    )

Function that maps the `f` function into all pairs of particles of
`system` that are found to be within the `cutoff`. 

The function `f` must be of the general form:
```
function f(x,y,i,j,d2,output)
    # operate on particle coordinates, distance and indexes
    # update output
    return output
end
```
where `x` and `y` are the coordinates (adjusted for the minimum
image) of the two particles involved, `i` and `j` their indices in the
original arrays of positions, `d2` their squared Euclidean distance,
and `output` the current value of the `output` variable. The `output`
variable must be updated within this function with the contribution
of the two particles involved. 

Thread-safety is taken care automatically in parallel executions.

`map_pairwise` is an alias to `map_pairwise!` for syntax consistency
when the `output` variable is immutable.

If `update_lists` is `false`, the cell lists will not be recomputed,
this may be useful for computing a different function from the same
coordinates.

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

julia> map_pairwise((x,y,i,j,d2,output) -> output += 1 / (1 + sqrt(d2)), sys)
1870.0274887950268
```

"""
function map_pairwise!(
    f::F,
    sys::AbstractParticleSystem;
    update_lists::Bool=true,
    show_progress::Bool=false
) where {F<:Function}
    sys.output = _reset_all_output!(sys.output, sys._output_threaded)
    UpdateParticleSystem!(sys, update_lists)
    sys.output = CellListMap.map_pairwise!(
        f, sys.output, sys._box, sys._cell_list;
        output_threaded=sys._output_threaded,
        parallel=sys.parallel,
        reduce=(output, output_threaded) -> reduce_output!(reducer, output, output_threaded),
        show_progress=show_progress
    )
    return sys.output
end
