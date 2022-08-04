
struct PeriodicSystem{B,C,A,O}
    box::B
    cell_list::C
    aux::A
    output_threaded::O
    parallel::Bool
end

copy_output(x::Number) = copy(x)
copy_output(x::Vector{<:Number}) = copy(x)
copy_output(x::Vector{<:SVector}) = copy(x)
copy_output(x::Vector{<:Tuple}) = copy(x)

function PeriodicSystem(
    x::AbstractVecOrMat, 
    sides::AbstractVector, 
    cutoff::Real,
    output;
    parallel::Bool = true,
    nbatches::Tuple{Int,Int} = (0,0),
)
    box = Box(sides, cutoff)
    cell_list = CellList(x, box; parallel = parallel, nbatches=nbatches)
    aux = AuxThreaded(cell_list)
    output_threaded = [ copy_output(output) for _ in 1:CellListMap.nbatches(cell_list) ] 
    sys = PeriodicSystem(box, cell_list, aux, output_threaded, parallel)
    return sys
end

import Base.show
function Base.show(io::IO,mime::MIME"text/plain",sys::PeriodicSystem)
    N = size(sys.box.unit_cell.matrix,1)
    println(io,"PeriodicSystem, containting:")
    show(io, mime, sys.box)
    println()
    show(io, mime, sys.cell_list)
    println("    Parallelization auxiliary data set for: ")
    show(io, mime, sys.cell_list.nbatches)
end