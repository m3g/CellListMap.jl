using CellListMap
using StaticArrays
using LinearAlgebra: norm
import Random

#
# In this example we test the implementation in a simple packmol-like 
# soft potential
#
function packmol(;parallel=false,sides=[46.4,31.5,18.7],tol=2.0,UnitCellType=OrthorhombicCell)

    function wrap(x,side)
        x = rem(x,side)
        if x >= side/2
            x -= side
        elseif x < -side/2
            x += side
        end
        return x
    end

    function mind(x,sides,tol)
        d = +Inf
        u = 0.
        for i in 1:length(x)-1
            for j in i+1:length(x)
                dd = norm(wrap.(x[i]-x[j],sides))
                d = min(d,dd)
                if dd < tol
                    u += (dd-tol)^2
                end
            end
        end
        return d, u
    end

    function fpair_cl(x,y,i,j,d2,f,box::Box)
        Δv = y - x
        d = sqrt(d2)
        fₓ = 2*(d - box.cutoff)*(Δv/d)
        f[i] += fₓ
        f[j] -= fₓ
        return f
    end

    function forces_cl!(f::Vector{T},x,box::Box,cl::CellList,fpair::F) where {T,F}
        fill!(f,zero(T))
        cl = UpdateCellList!(x,box,cl,parallel=parallel)
        map_pairwise!(
            (pair,f) -> fpair(pair.x,pair.y,pair.i,pair.j,pair.d2,f,box),
            f, box, cl, parallel=parallel
        )
        return f
    end

    function u_pack(x,box::Box,cl::CellList)
        cl = UpdateCellList!(x,box,cl,parallel=false)
        u = map_pairwise(
            (pair,u) -> begin
                u += (sqrt(pair.d2) - box.cutoff)^2 # objective function
                return u
            end,
            0., box, cl,
            parallel=parallel
        )
        return u
    end

    Random.seed!(321)

    x = [ sides .* rand(SVector{3,Float64}) for _ in 1:5000 ]

    mind_check, u_check = mind(x,sides,tol)

    box = Box(sides,tol,UnitCellType=UnitCellType)
    cl = CellList(x,box,parallel=parallel)
    u_pack = u_pack(x,box,cl)

    return u_pack ≈ u_check, u_pack, u_check, mind_check

end
