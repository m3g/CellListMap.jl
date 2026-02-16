using CellListMap
using StaticArrays
using LinearAlgebra: norm
import Random

#
# In this example we test the implementation in a simple packmol-like soft potential
#
function packmol(; parallel = false, sides = [46.4, 31.5, 18.7], tol = 2.0, UnitCellType = CellListMap.OrthorhombicCell)

    function wrap(x, side)
        x = rem(x, side)
        if x >= side / 2
            x -= side
        elseif x < -side / 2
            x += side
        end
        return x
    end

    function mind(x, sides, tol)
        d = +Inf
        u = 0.0
        for i in 1:(length(x) - 1)
            for j in (i + 1):length(x)
                dd = norm(wrap.(x[i] - x[j], sides))
                d = min(d, dd)
                if dd < tol
                    u += (dd - tol)^2
                end
            end
        end
        return d, u
    end

    Random.seed!(321)

    x = [ sides .* rand(SVector{3, Float64}) for _ in 1:5000 ]

    mind_check, u_check = mind(x, sides, tol)

    sys = ParticleSystem(
        positions = x,
        unitcell = sides,
        cutoff = tol,
        parallel = parallel,
        output = 0.0,
    )
    u = pairwise!(
        (pair, u) -> begin
            u += (pair.d - tol)^2
            return u
        end,
        sys
    )

    return u ≈ u_check, u, u_check, mind_check

end
