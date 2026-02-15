#
# Wrapper of the list of neighbors, that allows in-place updating of the lists
#
mutable struct NeighborList{T}
    n::Int
    list::Vector{Tuple{Int, Int, T}}
end

# Functions for parallel construction and reduction
function reset_output!(nb::NeighborList)
    nb.n = 0
    return nb
end
copy_output(nb::NeighborList) = NeighborList(nb.n, copy(nb.list))
function reducer!(nb1::NeighborList, nb2::NeighborList)
    ntot = nb1.n + nb2.n
    if length(nb1.list) < ntot
        resize!(nb1.list, ntot)
    end
    nb1.list[nb1.n + 1:ntot] .= nb2.list[1:nb2.n]
    nb1.n = ntot
    return nb1
end

# Function adds pair to the list
function push_pair!(pair, nb::NeighborList)
    (; i, j, d) =  pair
    nb.n += 1
    if nb.n > length(nb.list)
        push!(nb.list, (i, j, d))
    else
        nb.list[nb.n] = (i, j, d)
    end
    return nb
end

