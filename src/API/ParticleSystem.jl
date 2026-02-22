# Conversion helpers: Val{N} first so dispatch is on x/uc type alone, avoiding
# ambiguities with the SVector identity overloads.
_to_svectors(::Val{N}, x::AbstractVector{SVector{N,T}}) where {N,T} = x
_to_svectors(::Val{N}, x::AbstractVector{<:AbstractVector}) where {N} =
    Vector{SVector{N,eltype(eltype(x))}}([SVector{N,eltype(eltype(x))}(ntuple(j -> x[i][j], Val(N))) for i in eachindex(x)])
_to_svectors(::Val{N}, x::AbstractMatrix{T}) where {T,N} =
    Vector{SVector{N,T}}([SVector{N,T}(@view(x[:,i])) for i in axes(x, 2)])

_to_static_unitcell(::Val{N}, uc::SVector{N,T}) where {N,T} = uc
_to_static_unitcell(::Val{N}, uc::SMatrix{N,N,T,M}) where {N,T,M} = uc
_to_static_unitcell(::Val, ::Nothing) = nothing
_to_static_unitcell(::Val{N}, uc::AbstractVector) where {N} = SVector{N,eltype(uc)}(uc)
_to_static_unitcell(::Val{N}, uc::AbstractMatrix) where {N} = SMatrix{N,N,eltype(uc),N*N}(uc)

"""
    ParticleSystem(
        xpositions::AbstractVector{SVector{N,T}};
        ypositions::Union{AbstractVector{SVector{N,T}}, Nothing} = nothing,
        unitcell::Union{SVector, SMatrix, Nothing} = nothing,
        cutoff::Number,
        output,
        output_name::Symbol = :output,
        parallel::Bool = true,
        nbatches::Tuple{Int,Int} = (0, 0),
        lcell = 1,
        validate_coordinates = _validate_coordinates,
    )

Type-stable constructor. Accepts a `Vector{SVector{N,T}}` for positions and
an `SVector` (orthorhombic), `SMatrix` (triclinic), or `nothing` (non-periodic)
for the unit cell. All type parameters are resolved at compile time.

    ParticleSystem(;
        xpositions::Union{AbstractVector{<:AbstractVector},AbstractMatrix},
        #or
        xpositions::Union{AbstractVector{<:AbstractVector},AbstractMatrix},
        ypositions::Union{AbstractVector{<:AbstractVector},AbstractMatrix},
        # and
        unitcell::Union{Nothing,AbstractVecOrMat} = nothing,
        cutoff::Number,
        output,
        output_name::Symbol = :output,
        parallel::Bool = true,
        nbatches::Tuple{Int,Int} = (0, 0),
        validate_coordinates = _validate_coordinates,
    )

Flexible keyword-only constructor. Accepts any supported coordinate type
(vectors of vectors, matrices, or `SVector` arrays) and any `AbstractVecOrMat`
for the unit cell. Converts inputs and delegates to the type-stable constructor.

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

After construction, use `update!(system; xpositions=..., ypositions=..., cutoff=...,
unitcell=..., parallel=...)` to update any system properties before subsequent
`pairwise!` calls.

The `validate_coordinates` function can be used to validate the coordinates
before computations, and throw appropriate error messages. By default the validation checks if
the coordinates are not missing or NaN. The function must have a single
input parameter and be `(x) -> nothing` to skip any validation.

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
function ParticleSystem(
        xpositions::AbstractVector{SVector{N,T}};
        ypositions::Union{AbstractVector{SVector{N,T}}, Nothing} = nothing,
        unitcell::Union{SVector, SMatrix, Nothing} = nothing,
        cutoff::Number,
        output,
        output_name::Symbol = :output,
        parallel::Bool = true,
        nbatches::Tuple{Int,Int} = (0, 0),
        lcell = 1,
        validate_coordinates::F = _validate_coordinates,
    ) where {N, T, F<:Function}
    if isnothing(ypositions)
        _x = ParticleSystemPositions(xpositions)
        _uc = isnothing(unitcell) ? limits(_x; validate_coordinates) : unitcell
        _box = Box(_uc, cutoff; lcell)
        _cell_list = CellList(_x, _box; parallel, nbatches, validate_coordinates)
        _aux = AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list, :map)]
        output = _reset_all_output!(output, _output_threaded; reset = false)
        return ParticleSystem1{output_name}(_x, output, _box, _cell_list, _output_threaded, _aux, parallel, validate_coordinates)
    else
        _x = ParticleSystemPositions(xpositions)
        _y = ParticleSystemPositions(ypositions)
        _uc = isnothing(unitcell) ? limits(_x, _y; validate_coordinates) : unitcell
        _box = Box(_uc, cutoff; lcell)
        _cell_list = CellList(_x, _y, _box; parallel, nbatches, validate_coordinates)
        _aux = AuxThreaded(_cell_list)
        _output_threaded = [copy_output(output) for _ in 1:CellListMap.nbatches(_cell_list, :map)]
        output = _reset_all_output!(output, _output_threaded; reset = false)
        return ParticleSystem2{output_name}(_x, _y, output, _box, _cell_list, _output_threaded, _aux, parallel, validate_coordinates)
    end
end

# Thin wrapper: accepts any supported coordinate/unitcell types, converts, and
# delegates to the type-stable positional constructor above.
function ParticleSystem(;
        positions::SupportedCoordinatesTypes = nothing,
        xpositions::SupportedCoordinatesTypes = nothing,
        ypositions::SupportedCoordinatesTypes = nothing,
        unitcell::Union{AbstractVecOrMat, Nothing} = nothing,
        cutoff::Number,
        output,
        output_name::Symbol = :default_output_name,
        parallel::Bool = true,
        nbatches::Tuple{Int, Int} = (0, 0),
        lcell = 1,
        validate_coordinates::F = _validate_coordinates,
    ) where {F<:Function}
    if (isnothing(positions) && isnothing(xpositions)) || (!isnothing(positions) && !isnothing(xpositions))
        throw(ArgumentError("Either `positions` OR `xpositions` must be defined."))
    end
    xpositions = isnothing(positions) ? xpositions : positions
    DIM = get_dim(unitcell, xpositions, ypositions)
    vdim = Val(DIM)
    x  = _to_svectors(vdim, xpositions)
    y  = isnothing(ypositions) ? nothing : _to_svectors(vdim, ypositions)
    uc = _to_static_unitcell(vdim, unitcell)
    return ParticleSystem(x; ypositions=y, unitcell=uc,
        cutoff, output, output_name, parallel, nbatches, lcell, validate_coordinates)
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
setproperty!(sys::AbstractParticleSystem, ::Val{:unitcell}, x) = _update_unitcell!(sys, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:cutoff}, x) = _update_cutoff!(sys, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:parallel}, x) = setfield!(sys, :parallel, x)
# private properties
setproperty!(sys::AbstractParticleSystem, ::Val{:_box}, x) = setfield!(sys, :_box, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:_cell_list}, x) = setfield!(sys, :_cell_list, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:output}, x) = setfield!(sys, :output, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:_aux}, x) = setfield!(sys, :_aux, x)
setproperty!(sys::AbstractParticleSystem, ::Val{:_output_threaded}, x) = setfield!(sys, :_output_threaded, x)
