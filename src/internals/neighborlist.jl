#
# Wrapper of the list of neighbors, that allows in-place updating of the lists
#
mutable struct NeighborList{T}
    n::Int
    list::Vector{Tuple{Int, Int, T}}
end

import Base: push!, empty!, resize!, copy
empty!(x::NeighborList) = x.n = 0
function push!(x::NeighborList, pair)
    x.n += 1
    if x.n > length(x.list)
        push!(x.list, pair)
    else
        x.list[x.n] = pair
    end
    return x
end
function resize!(x::NeighborList, n::Int)
    x.n = n
    resize!(x.list, n)
    return x
end
copy(x::NeighborList{T}) where {T} = NeighborList{T}(x.n, copy(x.list))

# Function adds pair to the list
function push_pair!(i, j, d2, list::NeighborList)
    d = sqrt(d2)
    push!(list, (i, j, d))
    return list
end

# We have to define our own reduce function here (for the parallel version)
# this reduction can be dum assynchronously on a preallocated array
function reduce_lists(list::NeighborList{T}, list_threaded::Vector{<:NeighborList{T}}) where {T}
    ranges = cumsum(nb.n for nb in list_threaded)
    npairs = ranges[end]
    # need to resize here for the case where length(list) < npairs
    list = resize!(list, npairs)
    @sync for it in eachindex(list_threaded)
        lt = list_threaded[it]
        range = (ranges[it] - lt.n + 1):ranges[it]
        @spawn list.list[range] .= @view(lt.list[1:lt.n])
    end
    return list
end
