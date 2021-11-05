using CellListMap
using StaticArrays
using NearestNeighbors
using Test

function nl_NN(x,y,r)
    balltree = BallTree(x)
    return inrange(balltree,y,r,true)
end

function compare_nb_lists(list_CL,list_NN;self=false)
    if self 
        for (i,list_original) in pairs(list_NN)
            list = filter(!isequal(i),list_original)
            cl = filter(tup -> (tup[1] == i || tup[2] == i), list_CL)
            if length(cl) != length(list) 
                @show i
                @show length(list), list
                @show length(cl), cl
                return false
            end
            for j in list
                length(findall(tup -> (tup[1] == j || tup[2] == j), cl)) == 1 || return false
            end
        end
    else
        for (i,list) in pairs(list_NN)
            cl = filter(tup -> tup[2] == i, list_CL)
            if length(cl) != length(list) 
                @show i
                @show length(list), list
                @show length(cl), cl
                return false
            end
            for j in list
                length(findall(tup -> tup[1] == j, cl)) == 1 || return false
            end
        end
    end
    return true
end

@testset "neighbourlist" begin

    r = 0.1

    for N in [2,3]

        #
        # Using vectors as input
        #

        # With y smaller than x
        x = [ rand(SVector{N,Float64}) for _ in 1:1000 ]
        y = [ rand(SVector{N,Float64}) for _ in 1:500 ]

        nb = nl_NN(x,x,r)
        cl = CellListMap.neighbourlist(x,r)
        @test compare_nb_lists(cl,nb,self=true)

        nb = nl_NN(x,y,r)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=false)
        @test compare_nb_lists(cl,nb,self=false)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=true)
        @test compare_nb_lists(cl,nb,self=false)

        # with x smaller than y
        x = [ rand(SVector{N,Float64}) for _ in 1:500 ]
        y = [ rand(SVector{N,Float64}) for _ in 1:1000 ]
        nb = nl_NN(x,y,r)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=false)
        @test compare_nb_lists(cl,nb,self=false)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=true)
        @test compare_nb_lists(cl,nb,self=false)

        # Using matrices as input
        x = rand(N,1000)
        y = rand(N,500)

        nb = nl_NN(x,x,r)
        cl = CellListMap.neighbourlist(x,r)
        @test compare_nb_lists(cl,nb,self=true)

        nb = nl_NN(x,y,r)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=false)
        @test compare_nb_lists(cl,nb,self=false)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=true)
        @test compare_nb_lists(cl,nb,self=false)

        # with x smaller than y
        x = rand(N,500)
        y = rand(N,1000)
        nb = nl_NN(x,y,r)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=false)
        @test compare_nb_lists(cl,nb,self=false)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=true)
        @test compare_nb_lists(cl,nb,self=false)

        # Check random coordinates to test the limits more thoroughly
        check_random_NN = true
        for i in 1:500
            x = rand(SVector{N,Float64},100); y = rand(SVector{N,Float64},50); 
            nb = nl_NN(x,y,r); cl = CellListMap.neighbourlist(x,y,r,autoswap=false);
            check_random_NN = compare_nb_lists(cl,nb,self=false)
        end
        @test check_random_NN

        # with different types
        x = rand(Float32,N,500)
        y = rand(Float32,N,1000)
        nb = nl_NN(x,y,r)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=false)
        @test compare_nb_lists(cl,nb,self=false)
        cl = CellListMap.neighbourlist(x,y,r,autoswap=true)
        @test compare_nb_lists(cl,nb,self=false)

    end

end

