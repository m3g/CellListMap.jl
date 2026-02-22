#
# some auxiliary functions for testing neighbor lists
#
@testmodule TestingNeighborLists begin

    using CellListMap
    using LinearAlgebra: norm
    using StaticArrays

    export nl_NN
    export compare_nb_lists
    export is_unique
    export x_rotation, y_rotation, z_rotation, random_rotation

    function nl_NN(BallTree, inrange, x, y, r)
        balltree = BallTree(y)
        return inrange(balltree, x, r, true)
    end

    _norm(x::AbstractVector, y::AbstractVector, i, j) = norm(x[i] - y[j])
    _norm(x::AbstractMatrix, y::AbstractMatrix, i, j) = norm(x[:, i] - y[:, j])

    # for nb lists in a single set
    function compare_nb_lists(list_CL, list_NN, x, r::AbstractFloat)
        for (i, j_list) in pairs(list_NN)
            for j in j_list
                if i == j # inrange will return self-pairs
                    continue
                end
                ij_pairs = findall(p -> ((i, j) == (p[1], p[2])) || ((i, j) == (p[2], p[1])), list_CL)
                if length(ij_pairs) > 1
                    println("Non-unique pair: ", join((i, j, ij_pairs), " "))
                    return false, x, r
                end
                if length(ij_pairs) == 0
                    if !isapprox(_norm(x, x, i, j), r; atol = 1.0e-5)
                        println("Pair not found: ", join((i, j, ij_pairs, _norm(x, x, i, j)), " "))
                        return false, x, r
                    else
                        println("Warning: pair not found with d = r: ", join((i, j, _norm(x, x, i, j)), " "))
                    end
                end
            end
        end
        return true, nothing, nothing
    end


    # for nb lists in a single set
    function compare_nb_lists(list_CL, list_NN, x, y, r::AbstractFloat)
        for (i, j_list) in pairs(list_NN)
            for j in j_list
                ij_pairs = findall(p -> (i, j) == (p[1], p[2]), list_CL)
                if length(ij_pairs) > 1
                    println("Non-unique pair: ", join((i, j, j_list, list_CL[ij_pairs]), " "))
                    return false, x, y, r
                end
                if length(ij_pairs) == 0
                    if !isapprox(_norm(x, y, i, j), r; atol = 1.0e-5)
                        println("Pair not found: ", join((i, j, ij_pairs, _norm(x, y, i, j)), " "))
                        return false, x, y, r
                    else
                        println("Warning: pair not found with d = r: ", join((i, j, _norm(x, y, i, j)), " "))
                    end
                end
            end
        end
        return true, nothing, nothing, nothing
    end

    is_unique(list; self::Bool) = self ? is_unique_self(list) : is_unique_cross(list)
    is_unique_cross(list) = length(list) == length(unique(p -> (p[1], p[2]), list))
    is_unique_self(list) = length(list) == length(unique(p -> p[1] < p[2] ? (p[1], p[2]) : (p[2], p[1]), list))

    # Functions that define rotations along each axis, given the angle in 3D
    x_rotation(x) = @SMatrix[1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)]
    y_rotation(x) = @SMatrix[cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)]
    z_rotation(x) = @SMatrix[cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1]
    random_rotation() = z_rotation(2π * rand()) * y_rotation(2π * rand()) * x_rotation(2π * rand())

    # Test function to check if a pair in one list is unique in other list
    function _list_match_unique(lists_match, pair1, index_pair2, list2, cutoff; io, atol, first = "first", second = "second")
        # check if pair is unique
        if length(index_pair2) > 1
            println(io, "Non-unique pair: ", join((pair1, list2[index_pair2]), " "))
            println(io, "    On $first list: ", pair1)
            println(io, "    On $second list: ", list2[index_pair2])
            lists_match = false
        end
        # check if missing pair is at the cutof distance
        if length(index_pair2) == 0
            if !(isapprox(pair1[3] - cutoff, 0, atol = atol))
                println(io, "Pair of $first list not found in $second list: ", pair1)
                lists_match = false
            else
                println(io, "Warning: pair of $first list not found in $second list with d ≈ r: ", pair1)
            end
        end
        return lists_match
    end

    #=
    function lists_match(
        list1::Vector{Tuple{Int,Int,T1}},
        list2::Vector{Tuple{Int,Int,T2}},
        cutoff::Real;
        verbose::Bool=false,
        atol::Real=1e-10
    ) where {T1 <: Real, T2 <: Real}

# Extended help

Check if two neighbor lists are the same, up to a tolerance, and correctly taking into
account the possibility of pairs appearing or not if they are at the cutoff distance

## Example

```julia-repl
julia> using CellListMap

julia> x = [ rand(3) for i in 1:10_000 ];

julia> list1 = neighborlist(positions=x, cutoff=0.05);

julia> list2 = copy(list1)

julia> CellListMap.TestingNeighborLists.lists_match(list1, list2, 0.05; verbose = true)
true

julia> push!(list2, (1, 2, 0.05));

julia> CellListMap.TestingNeighborLists.lists_match(list1, list2, 0.05; verbose = true)
Warning: pair of second list not found in first list with d ≈ r: (1, 2, 0.05)
true

julia> push!(list2, (1, 3, 0.04));

julia> CellListMap.TestingNeighborLists.lists_match(list1, list2, 0.05; verbose = true)
Warning: pair of second list not found in first list with d ≈ r: (1, 2, 0.05)
Pair of second list not found in first list: (1, 3, 0.04)
false
```

=#
    function lists_match(
            list1::Vector{Tuple{Int, Int, T1}},
            list2::Vector{Tuple{Int, Int, T2}},
            cutoff::Real;
            verbose::Bool = false,
            atol::Real = 1.0e-10
        ) where {T1 <: Real, T2 <: Real}
        io = verbose ? stdout : devnull
        lists_match = true
        for pair1 in list1
            index_pair2 = findall(
                p -> ((p[1], p[2]) == (pair1[1], pair1[2])) | ((p[1], p[2]) == (pair1[2], pair1[1])), list2
            )
            # Check for uniqueness of pairs
            lists_match = _list_match_unique(lists_match, pair1, index_pair2, list2, cutoff; io, atol)
            # If pair is unique, check if distances match
            if length(index_pair2) == 1
                pair2 = list2[index_pair2[1]]
                # check if distances match
                if !isapprox(pair1[3] - pair2[3], 0, atol = atol)
                    println(io, "Warning: distances do not match: ", pair1, pair2)
                    lists_match = false
                end
            end
        end
        for pair2 in list2
            index_pair1 = findall(
                p -> ((p[1], p[2]) == (pair2[1], pair2[2])) | ((p[1], p[2]) == (pair2[2], pair2[1])), list1
            )
            # Check for uniqueness of pairs
            lists_match = _list_match_unique(lists_match, pair2, index_pair1, list1, cutoff; io, atol, first = "second", second = "first")
        end
        return lists_match
    end

end # module
