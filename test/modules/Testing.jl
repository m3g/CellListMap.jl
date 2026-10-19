@testmodule Testing begin

    using CellListMap
    using StaticArrays
    using Setfield: @set!
    using CellListMap: ParticleSystem1, ParticleSystem2
    using CellListMap: reset_output!

    export pathological_coordinates
    export map_naive!
    export check_random_cells
    export test_random_cells
    export xatomic, xgalactic
    export test_pathological
    export missing_pair
    export simple_test
    export PSP
    export f1, f2, _f1, _f2

    const PSP = CellListMap.ParticleSystemPositions

    f1(pair, out) = out += pair.d 
    f2(pair, out) = out += pair.d2 
    _f1(x, y, i, j, d2, out) = out += sqrt(d2)
    _f2(x, y, i, j, d2, out) = out += d2

    #=
    pathological_coordinates(N)

# Extended help

Function to generate some coordinates with pathological properties, for testing.
Returns `x`, `y`, `sides` and `cutoff`.

=#
    function pathological_coordinates(N)

        # Number of particles, sides and cutoff
        sides = @SVector [250.0, 250.0, 250.0]
        cutoff = 10.0

        # Particle positions
        x = [sides .* rand(SVector{3, Float64}) for i in 1:N]
        # Add some pathological coordinates
        x[1] = -sides / 2
        x[2] = -sides / 2 + @SVector [nextfloat(0.0), nextfloat(0.0), sides[3] * rand()]
        x[3] = -sides / 2 + @SVector [prevfloat(0.0), prevfloat(0.0), sides[3] * rand()]
        x[4] = sides / 2 + @SVector [nextfloat(0.0), nextfloat(0.0), sides[3] * rand()]
        x[5] = sides / 2 + @SVector [prevfloat(0.0), prevfloat(0.0), sides[3] * rand()]
        x[10] = sides
        x[11] = sides + @SVector [nextfloat(0.0), nextfloat(0.0), sides[3] * rand()]
        x[12] = sides + @SVector [prevfloat(0.0), prevfloat(0.0), sides[3] * rand()]
        x[13] = @SVector [nextfloat(0.0), nextfloat(0.0), sides[3] * rand()]
        x[14] = @SVector [prevfloat(0.0), prevfloat(0.0), sides[3] * rand()]
        x[15] = zero(SVector{3, Float64})
        x[100] = @SVector [sides[1] / 2, -sides[2] / 2, 2 * sides[3]]
        y = [sides .* rand(SVector{3, Float64}) for i in 1:N]

        return x, y, sides, cutoff
    end

    map_naive!(f, sys::ParticleSystem1) =
        map_naive!(f, reset_output!(sys.output), sys.xpositions.x, sys.private.box)
    map_naive!(f, sys::ParticleSystem2) =
        map_naive!(f, reset_output!(sys.output), sys.xpositions.x, sys.ypositions.x, sys.private.box)

    #=
    map_naive!(f::Function, output, x::AbstractVector, box::Box)
    map_naive!(f::Function, output, x::AbstractVector, y::AbstractVector, box::Box)

    # Extended help
    
    Function that uses the naive pairwise mapping algorithm, for testing.
    
    =#
    function map_naive!(f::Function, output, x::AbstractVector, box::CellListMap.Box)
        (; input_unit_cell, cutoff_sqr) = box
        for i in 1:(length(x) - 1)
            xᵢ = x[i]
            for j in (i + 1):length(x)
                xⱼ = CellListMap.wrap_relative_to(x[j], xᵢ, input_unit_cell.matrix)
                d2 = sum(abs2, xᵢ - xⱼ)
                if d2 <= cutoff_sqr
                    output = f(xᵢ, xⱼ, i, j, d2, output)
                end
            end
        end
        return output
    end

    function map_naive!(f::Function, output, x::AbstractVector, y::AbstractVector, box::CellListMap.Box)
        (; input_unit_cell, cutoff_sqr) = box
        for i in eachindex(x)
            xᵢ = x[i]
            for j in eachindex(y)
                yⱼ = CellListMap.wrap_relative_to(y[j], xᵢ, input_unit_cell.matrix)
                d2 = sum(abs2, xᵢ - yⱼ)
                if d2 <= cutoff_sqr
                    output = f(xᵢ, yⱼ, i, j, d2, output)
                end
            end
        end
        return output
    end


    test_map2(box, cl; parallel = true) =
        CellListMap._pairwise!(
        (pair, s) -> begin
            s += 1
        end,
        0, box, cl, parallel = parallel
    )

    function missing_pair(x1, x2)
        for x in eachrow(x2)
            found = false
            for y in eachrow(x1)
                (x[1] == y[1]) && (x[2] == y[2]) && (found = true)
                (x[1] == y[2]) && (x[2] == y[1]) && (found = true)
                if found
                    break
                end
            end
            if !found
                return x
            end
        end
    end

    test_map(box, cl; parallel = true) = CellListMap._pairwise!((pair, s) -> s += 1 / pair.d2, 0.0, box, cl, parallel = parallel)
    test_naive(x, box) = map_naive!((x, y, i, j, d2, s) -> s += 1 / d2, 0.0, x, box)
    function simple_test(x, box; parallel = true)
        cl = CellListMap.CellList(PSP(x), box, parallel = parallel)
        r_naive = test_naive(x, box)
        r_map = test_map(box, cl, parallel = parallel)
        println("naive = $r_naive")
        println("map   = $r_map")
        return "Test passed: $(r_naive ≈ r_map)"
    end

    function check_random_cells(
            N = 3, M = 10;
            show_progress = true,
            lcell = 1,
            parallel = false,
            UnitCellType = CellListMap.OrthorhombicCell
        )
        local x, box
        ntrial = 0
        show_progress && (p = Progress(10000, 0.5))
        while ntrial < 2000
            unit_cell_matrix = zeros(SMatrix{N, N, Float64})
            if UnitCellType == CellListMap.OrthorhombicCell
                for i in 1:N
                    @set! unit_cell_matrix[i, i] = -10 + 20 * rand()
                end
            else
                for i in 1:N, j in 1:N
                    @set! unit_cell_matrix[i, j] = -10 + 20 * rand()
                end
            end
            cutoff = 1 + rand()
            try
                if !(CellListMap.check_unit_cell(unit_cell_matrix, cutoff))
                    continue
                end
            catch
                return false, unit_cell_matrix, cutoff
            end
            box = try
                if UnitCellType == CellListMap.OrthorhombicCell
                    box = CellListMap.Box([unit_cell_matrix[i, i] for i in 1:N], cutoff, lcell = lcell)
                else
                    box = CellListMap.Box(unit_cell_matrix, cutoff, lcell = lcell)
                end
            catch
                return false, CellListMap.UnitCellType, (unit_cell_matrix, cutoff, lcell)
            end
            if prod(box.nc) > 100000
                continue
            end
            show_progress && (next!(p))
            x = 10 * rand(SVector{N, Float64}, M)
            for i in eachindex(x)
                x[i] = x[i] .- 50
            end
            test = try
                cl = CellListMap.CellList(PSP(x), box, parallel = parallel)
                test_map(box, cl, parallel = parallel)
            catch
                return false, x, box
            end
            if test ≈ 0
                continue
            end
            ntrial += 1
            if !(test ≈ test_naive(x, box))
                show_progress && println("FOUND PROBLEMATIC SETUP.")
                return false, x, box
            else

            end
            if parallel
                aux = CellListMap.AuxThreaded(cl)
                cl = CellListMap.UpdateCellList!(PSP(x), box, cl, aux, parallel = parallel)
            else
                cl = CellListMap.UpdateCellList!(PSP(x), box, cl, parallel = parallel)
            end
            if !(test_map(box, cl, parallel = parallel) ≈ test_naive(x, box))
                show_progress && println("FOUND PROBLEMATIC SETUP.")
                return false, x, box
            end
        end
        show_progress && println(" ALL PASSED! ")
        return true, nothing, nothing
    end

    function test_random_cells(::Type{T}) where {T}
        for N in 2:3,
                M in rand(10:20),
                parallel in [false, true],
                lcell in 1:3

            test = check_random_cells(
                N, M,
                UnitCellType = T,
                parallel = parallel,
                lcell = lcell,
                show_progress = false
            )
            if test[1] != true
                return test
            end
        end
        return nothing, nothing, nothing
    end

    function drawbox(box::CellListMap.Box{UnitCellType, 2}) where {UnitCellType}
        S = SVector{2, Float64}
        m = box.input_unit_cell.matrix
        x = [
            S(0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 2]),
            S(m[:, 2]),
            S(0.0, 0.0),
        ]
        return x
    end

    function drawbox(box::CellListMap.Box{UnitCellType, 3}) where {UnitCellType}
        S = SVector{3, Float64}
        m = box.input_unit_cell.matrix
        x = [
            S(0.0, 0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 2]),
            S(m[:, 2]),
            S(0.0, 0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 3]),
            S(m[:, 3]),
            S(0.0, 0.0, 0.0),
            S(m[:, 2]),
            S(m[:, 2] + m[:, 3]),
            S(m[:, 2]),
            S(0.0, 0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 2]),
            S(m[:, 1] + m[:, 2]) + m[:, 3],
            S(m[:, 1] + m[:, 3]),
            S(m[:, 3]),
            S(m[:, 3]) + m[:, 2],
            S(m[:, 1] + m[:, 2]) + m[:, 3],
        ]
        return x
    end


    # Create a set of random particles with an atomic density similar to that
    # of water at room temperature and pressure.
    function xatomic(n = 100_000; cutoff = 12.0, lcell = 1, type = :vec)
        atomic_density_of_water = 0.1 # atoms/ Angstrom^3
        vol = n / atomic_density_of_water
        l = vol^(1 / 3)
        if type == :mat
            x = l * rand(3, n)
        elseif type == :vec
            x = [l * rand(SVector{3, Float64}) for _ in 1:n]
        end
        box = CellListMap.Box([l, l, l], cutoff, lcell = lcell)
        return x, box
    end

    # Create a set of  random particles with the typical density of galaxies
    function xgalactic(n = 8_000_000; cutoff = 5.0, lcell = 1, type = :vec)
        galactic_density = 12 # galaxyes per cubic mega-parsecs (https://doi.org/10.1093/mnras/106.2.121)
        vol = n / galactic_density
        l = vol^(1 / 3) # mega-parsecs
        if type == :mat
            x = l * rand(3, n)
        elseif type == :vec
            x = [l * rand(SVector{3, Float64}) for _ in 1:n]
        end
        box = CellListMap.Box([l, l, l], cutoff, lcell = lcell)
        return x, box
    end


    function test_pathological(; nmax = 1000, npoints = 2)

        function g(i, j, d2, cutoff, out)
            if !(d2 ≈ cutoff^2)
                out += i + j
            end
            return out
        end

        l = sqrt(2) / 2
        matrices = [
            @SMatrix[1 0; 0 1],
            @SMatrix[l 0; l 1],
            @SMatrix[1.1 0; 0 1],
            @SMatrix[1.2 0; 0 1],
            @SMatrix[1 0; 0 1.1],
            @SMatrix[1 0; 0 1.2],
            @SMatrix[1 0.2; 0 1.2],
            @SMatrix[1 0.2; 0.2 1.2],
            @SMatrix[1.2 0.2; 0.2 1.2],
            @SMatrix[-1.2 0.2; 0.2 1.2],
            @SMatrix[-1.2 0.2; 0.2 -1.2],
        ]

        for n in 0:nmax, matrix in matrices, lcell in 1:5
            x = rand(SVector{2, Float64}, npoints)
            cutoff = 0.2
            box = CellListMap.Box(matrix, cutoff; lcell = lcell)
            naive = map_naive!((x, y, i, j, d2, out) -> g(i, j, d2, box.cutoff, out), 0, x, box)
            cl = CellListMap.CellList(PSP(x), box)
            clmap = CellListMap._pairwise!((pair, out) -> g(pair.i, pair.j, pair.d2, box.cutoff, out), 0, box, cl)
            if naive != clmap
                return x, box
            end
        end
        return nothing, nothing
    end

    #
    # Functions for drawing unitcells, for testing
    #
    #=
    draw_cell_vertices(m::AbstractMatrix)

Function that returns the vertices of a unit cell matrix in 2D or 3D, as a vector
of static vectors, in a proper order for plotting the cell (the first vertex, in the
origin, is repeated at the end of the list, to close the figure)

=#
    function draw_cell_vertices(m::AbstractMatrix)
        if size(m) == (2, 2)
            x = _draw_cell_vertices2D(m)
        elseif size(m) == (3, 3)
            x = _draw_cell_vertices3D(m)
        else
            throw(ArgumentError("draw_cell_vertices only supports square matrices in 2 or 3 dimensions."))
        end
    end

    function _draw_cell_vertices2D(m::AbstractMatrix)
        S = SVector{2, Float64}
        x = [
            S(0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 2]),
            S(m[:, 2]),
            S(0.0, 0.0),
        ]
        return x
    end

    function _draw_cell_vertices3D(m::AbstractMatrix)
        S = SVector{3, Float64}
        x = [
            S(0.0, 0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 2]),
            S(m[:, 2]),
            S(0.0, 0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 3]),
            S(m[:, 3]),
            S(0.0, 0.0, 0.0),
            S(m[:, 2]),
            S(m[:, 2] + m[:, 3]),
            S(m[:, 2]),
            S(0.0, 0.0, 0.0),
            S(m[:, 1]),
            S(m[:, 1] + m[:, 2]),
            S(m[:, 1] + m[:, 2]) + m[:, 3],
            S(m[:, 1] + m[:, 3]),
            S(m[:, 3]),
            S(m[:, 3]) + m[:, 2],
            S(m[:, 1] + m[:, 2]) + m[:, 3],
        ]
        return x
    end

    #=
    draw_cell(m::AbstractMatrix; aspect_ratio=:auto)

Draw the unit cell in a 2D or 3D plot. Requires `using Plots`.

=#
    function draw_cell(m::AbstractMatrix; aspect_ratio = :auto)
        plot = Main.plot
        plot! = Main.plot!
        vertices = draw_cell_vertices(m)
        plt = plot(Tuple.(vertices), label = nothing)
        lims = cell_limits(m)
        dx = (lims[2][1] - lims[1][1]) / 10
        dy = (lims[2][2] - lims[1][2]) / 10
        plot!(
            plt,
            xlims = (lims[1][1] - dx, lims[2][1] + dx),
            ylims = (lims[1][2] - dy, lims[2][2] + dy),
            aspect_ratio = aspect_ratio,
            xlabel = "x",
            ylabel = "y",
            zlabel = "z",
        )
        return plt
    end

end # testing module
