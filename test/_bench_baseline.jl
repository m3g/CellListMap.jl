using CellListMap
using StaticArrays
using Chairmarks

f(pair, out) = out += pair.d

function xatomic(n)
    atomic_density_of_water = 0.1
    vol = n / atomic_density_of_water
    l = vol^(1 / 3)
    x = [l * rand(SVector{3, Float64}) for _ in 1:n]
    return x, l
end

function make_ortho_sys(n; cutoff = 12.0)
    x, l = xatomic(n)
    sys = ParticleSystem(
        positions = x,
        cutoff = cutoff,
        unitcell = [l, l, l],
        output = 0.0,
    )
    return sys
end

# Triclinic: skew the cell so it's clearly non-orthorhombic
function make_triclinic_sys(n; cutoff = 12.0)
    x, l = xatomic(n)
    # Upper-triangular cell matrix (valid triclinic)
    unitcell = [l 0.1*l 0.05*l; 0 l 0.1*l; 0 0 l]
    # Wrap coordinates into the triclinic cell
    x_wrapped = [mod.(xi, l) for xi in x]  # approximate wrap
    sys = ParticleSystem(
        positions = x_wrapped,
        cutoff = cutoff,
        unitcell = unitcell,
        output = 0.0,
    )
    return sys
end

# Cross-interaction benchmark
function make_cross_sys(n; cutoff = 12.0)
    x, l = xatomic(n)
    y = [l * rand(SVector{3, Float64}) for _ in 1:n÷2]
    sys = ParticleSystem(
        xpositions = x,
        ypositions = y,
        cutoff = cutoff,
        unitcell = [l, l, l],
        output = 0.0,
    )
    return sys
end

println("=== Baseline Benchmark ===")
println("Julia version: ", VERSION)
println()

for n in [10_000, 100_000]
    println("--- n = ", n, " ---")

    # Orthorhombic self
    sys = make_ortho_sys(n)
    t = @b pairwise!(f, sys)
    println("  ortho  self:  ", t.time * 1e3, " ms")

    # Triclinic self
    sys_tri = make_triclinic_sys(n)
    t = @b pairwise!(f, sys_tri)
    println("  tri    self:  ", t.time * 1e3, " ms")

    # Orthorhombic cross
    sys_cross = make_cross_sys(n)
    t = @b pairwise!(f, sys_cross)
    println("  ortho  cross: ", t.time * 1e3, " ms")

    println()
end
