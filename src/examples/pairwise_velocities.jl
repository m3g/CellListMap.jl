using CellListMap
using StaticArrays
import Random

#
# Here we simulate the computation of the histogram of pairwise velocities, for a set of 
# particles for which the distances are known. This histogram is useful in the study
# galaxy motions, for example. 
#
function pairwise_velocities(::Type{T}=Float64;
    N=100_000,cd=true,parallel=true,lcell=1,nbatches=zero(CellListMap.NumberOfBatches)
) where T

    @inline dot(x::SVector{3,T}, y::SVector{3,T}) where T = x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
    
    function compute_pairwise_mean_cell_lists!(x, y, i, j, d2, hist, velocities, rbins)
        d = x - y
        r = sqrt(d2)
        ibin = searchsortedfirst(rbins, r) - 1
        hist[1][ibin] += 1
        hist[2][ibin] += dot(velocities[i] - velocities[j], d) / r
        return hist
    end
  
    function reduce_hist(hist, hist_threaded)
        for i in 1:length(hist_threaded)
            hist[1] .+= hist_threaded[i][1]
            hist[2] .+= hist_threaded[i][2]
        end
        return hist
    end
  
    n_halos = N
  
    if cd
        density = T(10^5 / 250^3)  # density of the original problem
        boxsize = T((n_halos/density)^(1/3))
    else
        boxsize = T(250.)
    end

    Random.seed!(321)
    Lbox = T[boxsize,boxsize,boxsize]
    positions = convert.(T,boxsize .* rand(Float64, 3, n_halos))
    velocities = convert.(T,rand(Float64, 3, n_halos))
    rbins = T[0.,2.,4.,6.,8.,10.]
    r_max = maximum(rbins)
  
    n = size(positions)[2]
    positions = reshape(reinterpret(SVector{3,T}, positions), n)
    velocities = reshape(reinterpret(SVector{3,T}, velocities), n)
  
    box = Box(Lbox, r_max, UnitCellType=OrthorhombicCell, lcell=lcell) 
    cl = CellList(positions, box, parallel=parallel, nbatches=nbatches)
    hist = (zeros(Int, length(rbins) - 1), zeros(T, length(rbins) - 1))
  
    # Needs this to stabilize the type of velocities and hist, probably
    function barrier!(f!, velocities, rbins, hist, box, cl, reduce_hist, parallel)
        map_pairwise!(
            (x, y, i, j, d2, hist) -> f!(x, y, i, j, d2, hist, velocities, rbins),
            hist, box, cl,
            reduce=reduce_hist,
            parallel=parallel
        )
        return hist
    end
  
    barrier!(
        compute_pairwise_mean_cell_lists!,
        velocities,rbins,hist,box,cl,reduce_hist,parallel
    )
  
    n_pairs = hist[1]
    mean_v_r = hist[2]
    mean_v_r[n_pairs .> 0] = mean_v_r[n_pairs .> 0] ./ n_pairs[n_pairs .> 0]
  
    correct_result = [ 
         0.0007883474482652579
        -0.0035662371878635722
        -0.00040742008823982926
         0.0003623877989466509
        -0.0010441334538614498
    ]
    return mean_v_r ≈ correct_result, mean_v_r

end

# Credits to @florpi, and we use her github name here as an alias for testing
# and historical reasons.
const florpi = pairwise_velocities
