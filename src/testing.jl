"""
    pathological_coordinates(N)

$(INTERNAL)

# Extended help

Function to generate some coordinates with pathological properties, for testing.
Returns `x`, `y`, `sides` and `cutoff`.

"""
function pathological_coordinates(N)

    # Number of particles, sides and cutoff
    sides = @SVector [250.0, 250.0, 250.0]
    cutoff = 10.0

    # Particle positions
    x = [sides .* rand(SVector{3,Float64}) for i in 1:N]
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
    x[100] = @SVector [sides[1] / 2, -sides[2] / 2, 2 * sides[3]]
    y = [sides .* rand(SVector{3,Float64}) for i in 1:N]

    return x, y, sides, cutoff
end

"""
    map_naive!(f::Function, output, x::AbstractVector, box::Box)
    map_naive!(f::Function, output, x::AbstractVector, y::AbstractVector, box::Box)

$(INTERNAL)

# Extended help

Function that uses the naive pairwise mapping algorithm, for testing.

"""
function map_naive!(f::Function, output, x::AbstractVector, box::Box)
    @unpack input_unit_cell, cutoff_sqr = box
    for i in 1:length(x)-1
        xᵢ = x[i]
        for j in i+1:length(x)
            xⱼ = wrap_relative_to(x[j], xᵢ, input_unit_cell.matrix)
            d2 = norm_sqr(xᵢ - xⱼ)
            if d2 <= cutoff_sqr
                output = f(xᵢ, xⱼ, i, j, d2, output)
            end
        end
    end
    return output
end

function map_naive!(f::Function, output, x::AbstractVector, y::AbstractVector, box::Box)
    @unpack input_unit_cell, cutoff_sqr = box
    for i in eachindex(x)
        xᵢ = x[i]
        for j in eachindex(y)
            yⱼ = wrap_relative_to(y[j], xᵢ, input_unit_cell.matrix)
            d2 = norm_sqr(xᵢ - yⱼ)
            if d2 <= cutoff_sqr
                output = f(xᵢ, yⱼ, i, j, d2, output)
            end
        end
    end
    return output
end

"""
    view_celllist_particles(cl::CellList)

$(INTERNAL)

# Extended help

Auxiliary function to view the particles of a computing box, including images created
for computing purposes.

"""
function view_celllist_particles(cl::CellList{N,T}) where {N,T}
    x = SVector{N,T}[]
    for icell in 1:cl.n_cells_with_particles
        cell = cl.cells[icell]
        for ip in 1:cell.n_particles
            p = cell.particles[ip]
            push!(x, p.coordinates)
        end
    end
    return x
end

test_map2(box, cl; parallel=true) =
    map_pairwise!(
        (x, y, i, j, d2, s) -> begin
            s += 1
        end,
        0, box, cl, parallel=parallel
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

test_map(box, cl; parallel=true) = map_pairwise!((x, y, i, j, d2, s) -> s += 1 / d2, 0.0, box, cl, parallel=parallel)
test_naive(x, box) = CellListMap.map_naive!((x, y, i, j, d2, s) -> s += 1 / d2, 0.0, x, box)
function simple_test(x, box; parallel=true)
    cl = CellList(x, box, parallel=parallel)
    r_naive = test_naive(x, box)
    r_map = test_map(box, cl, parallel=parallel)
    println("naive = $r_naive")
    println("map   = $r_map")
    return "Test passed: $(r_naive ≈ r_map)"
end

function check_random_cells(
    N=3, M=10;
    show_progress=true,
    lcell=1,
    parallel=false,
    UnitCellType=OrthorhombicCell
)
    local x, box
    ntrial = 0
    show_progress && (p = Progress(10000, 0.5))
    while ntrial < 2000
        unit_cell_matrix = zeros(SMatrix{N,N,Float64})
        if UnitCellType == OrthorhombicCell
            for i in 1:N
                @set! unit_cell_matrix[i, i] = 1 + 10 * rand()
            end
        else
            for i in 1:N, j in 1:N
                @set! unit_cell_matrix[i, j] = 1 + 10 * rand()
            end
        end
        cutoff = 1 + rand()
        try
            if !check_unit_cell(unit_cell_matrix, cutoff, printerr=false)
                continue
            end
        catch
            return false, unit_cell_matrix, cutoff
        end
        box = try
            if UnitCellType == OrthorhombicCell
                box = Box([unit_cell_matrix[i, i] for i in 1:N], cutoff, lcell=lcell)
            else
                box = Box(unit_cell_matrix, cutoff, lcell=lcell)
            end
        catch
            return false, UnitCelType, (unit_cell_matrix, cutoff, lcell)
        end
        if prod(box.nc) > 100000
            continue
        end
        show_progress && (next!(p))
        x = 10 * rand(SVector{N,Float64}, M)
        for i in eachindex(x)
            x[i] = x[i] .- 50
        end
        test = try 
            cl = CellList(x, box, parallel=parallel)
            test_map(box, cl, parallel=parallel)
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
            aux = AuxThreaded(cl)
            cl = UpdateCellList!(x, box, cl, aux, parallel=parallel)
        else
            cl = UpdateCellList!(x, box, cl, parallel=parallel)
        end
        if !(test_map(box, cl, parallel=parallel) ≈ test_naive(x, box))
            show_progress && println("FOUND PROBLEMATIC SETUP.")
            return false, x, box
        end
    end
    show_progress && println(" ALL PASSED! ")
    return true, nothing, nothing
end

function test_random_cells()
    for N in 2:3, 
        M in rand(10:20), 
        UnitCellType in [ TriclinicCell, OrthorhombicCell ],
        parallel in [ false, true ],
        lcell in 1:3
        test = CellListMap.check_random_cells(
            N,M,
            UnitCellType=UnitCellType,
            parallel=parallel,
            lcell=lcell,
            show_progress=false
        )
        if test[1] != true
            return test
        end
    end
    return nothing, nothing, nothing
end

@testitem "random cells" begin
    using CellListMap
    using StaticArrays
    # Test random cells of all possible types
    @test CellListMap.test_random_cells() == (nothing, nothing, nothing) 
end

function drawbox(box::Box{UnitCellType,2}) where {UnitCellType}
    S = SVector{2,Float64}
    m = box.input_unit_cell.matrix
    x = [
        S(0.0, 0.0),
        S(m[:, 1]),
        S(m[:, 1] + m[:, 2]),
        S(m[:, 2]),
        S(0.0, 0.0)
    ]
    return x
end

function drawbox(box::Box{UnitCellType,3}) where {UnitCellType}
    S = SVector{3,Float64}
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
        S(m[:, 1] + m[:, 2]) + m[:, 3]
    ]
    return x
end

function get_particles(cl::CellList{N,T}) where {N,T}
    real = SVector{N,T}[]
    ghost = SVector{N,T}[]
    for cell in cl.cells
        for i in 1:cell.n_particles
            p = cell.particles[i]
            if p.real
                push!(real, p.coordinates)
            else
                push!(ghost, p.coordinates)
            end
        end
    end
    return real, ghost
end

"""
    draw_computing_cell(x,box::Box{UnitCellType,2}) where UnitCellType
    draw_computing_cell(cl::CellList,box::Box{UnitCellType,2},x) where UnitCellType

$(INTERNAL)

This function creates a plot of the computing cell, in two dimensions.

"""
function draw_computing_cell(x, box::Box{UnitCellType,2};
    parallel=true,
    xticks=nothing,
    yticks=nothing
) where {UnitCellType}
    cl = CellList(x, box, parallel=parallel)
    return draw_computing_cell(cl, box; xticks=xticks, yticks=yticks)
end
function draw_computing_cell(
    cl::CellList, box::Box{UnitCellType,2};
    xticks=nothing,
    yticks=nothing,
    x=nothing
) where {UnitCellType}
    real, ghost = get_particles(cl)
    plt = Main.plot()
    vertices = draw_cell_vertices(box.aligned_unit_cell.matrix)
    Main.plot!(plt, Tuple.(vertices), label=:none)
    Main.scatter!(plt, Tuple.(real), label=:none, color=:blue, markeralpha=1)
    Main.scatter!(plt, Tuple.(ghost), label=:none, color=:blue, markeralpha=0.3)
    if !isnothing(x)
        x_wrapped = wrap_to_first.(x, Ref(box.input_unit_cell.matrix))
        x_rotated = Ref(box.rotation) .* x_wrapped
        Main.scatter!(plt, Tuple.(wrap_to_first.(x_rotated, Ref(box.aligned_unit_cell.matrix))), color=:red, label=:none)
    end
    xmin, xmax = cell_limits(box.aligned_unit_cell.matrix)
    xmin = xmin .- 3*box.cell_size
    xmax = xmax .+ 3*box.cell_size
    isnothing(xticks) && (xticks = (round.(digits=3, xmin[1]:box.cell_size[1]:xmax[1])))
    isnothing(yticks) && (yticks = (round.(digits=3, xmin[2]:box.cell_size[2]:xmax[2])))
    Main.plot!(plt,
        aspect_ratio=1, framestyle=:box, xrotation=60,
        xlims=(xmin[1], xmax[1]),
        ylims=(xmin[2], xmax[2]),
        xticks=xticks,
        yticks=yticks,
        title="cell_sizes = $(box.cell_size)",
        titlefontsize=8,
    )
    return plt
end

"""
    draw_computing_cell(x,box::Box{UnitCellType,3}) where UnitCellType

$(INTERNAL)

This function creates a plot of the computing cell, in three dimensions.

"""
function draw_computing_cell(x, box::Box{UnitCellType,3}; parallel=true) where {UnitCellType}
    cl = CellList(x, box, parallel=parallel)
    vertices = draw_cell_vertices(box.aligned_unit_cell.matrix)
    p = view_celllist_particles(cl)
    plt = Main.plot()
    Main.plot!(plt, Tuple.(vertices), label=:none)
    Main.scatter!(plt, Tuple.(p), label=:none, markeralpha=0.3)
    x_rotated = Ref(box.rotation) .* x
    Main.scatter!(plt, Tuple.(wrap_to_first.(x_rotated, Ref(box.aligned_unit_cell.matrix))), label=:none, markeralpha=0.3)
    lims = Vector{Float64}[]
    xmin, xmax = box.computing_limits
    xmin = xmin .- box.cell_size
    xmax = xmax .+ box.cell_size
    Main.plot!(plt,
        aspect_ratio=1, framestyle=:box, xrotation=60, yrotation=-70, zrotation=0,
        xlims=(xmin[1], xmax[1]),
        ylims=(xmin[2], xmax[2]),
        zlims=(xmin[3], xmax[3]),
        xticks=(round.(digits=3, xmin[1]:box.cell_size[1]:xmax[1])),
        yticks=(round.(digits=3, xmin[2]:box.cell_size[2]:xmax[2])),
        zticks=(round.(digits=3, xmin[3]:box.cell_size[3]:xmax[3])),
    )
    return plt
end


function compare_cells(cl1::CellList{N,T}, cl2::CellList) where {N,T}

    if cl1.n_cells_with_real_particles != cl2.n_cells_with_real_particles
        println("Error: n_cells_with_real_particles differ")
    end
    if cl1.n_particles != cl2.n_particles
        println("Error: n_particles differ")
    end
    for i in 1:cl1.n_cells_with_particles
        linear_index =
            j = findfirst(
                isequal(cl1.cell_indices[i]),
                @view(cl2.cell_indices[1:cl2.n_cells_with_particles])
            )
        if isnothing(j)
            println(" Warning: cells_with_particles: Could not find index $(cl1.cell_indices[i])")
            break
        end
    end
    for i in 1:cl1.n_cells_with_real_particles
        j = findfirst(
            isequal(cl1.cell_indices_real[i]),
            @view(cl2.cell_indices_real[1:cl2.n_cells_with_real_particles])
        )
        if isnothing(j)
            println(" Warning: cells_with_real_particles: Could not find index $(cl1.cell_indices_real[i])")
            break
        end
    end
    if length(cl1.projected_particles) != length(cl2.projected_particles)
        println("Error: length of project_particles differ.")
    end
    ncheck = 0
    nparticle_check = 0
    if length(cl1.cells) != length(cl2.cells)
        println("Error: lengths of cells lists differ")
    else
        differ = false
        for i in 1:cl1.n_cells_with_particles
            if differ
                break
            end
            ncheck += 1
            ci = cl1.cells[i]
            j = findfirst(c -> c.linear_index == ci.linear_index, cl2.cells)
            if isnothing(j)
                println("Error: linear index of $(ci.linear_index) not found.")
                differ = true
            else
                cj = cl2.cells[j]
                if ci.cartesian_index != cj.cartesian_index
                    println("Error: cartesian_index differ")
                    differ = true
                end
                if !(ci.center ≈ cj.center)
                    println("Error: centers differ")
                    differ = true
                end
                if ci.contains_real != cj.contains_real
                    println("Error: contains_real differ")
                    differ = true
                end
                if ci.n_particles != cj.n_particles
                    println("Error: n_particles differ")
                    differ = true
                end
                for p1 in ci.particles
                    j = findfirst(p -> p.index == p1.index, cj.particles)
                    if isnothing(j)
                        println("Error: could not find particle in cl2.")
                        @show p1
                        differ = true
                    else
                        nparticle_check += 1
                        p2 = cj.particles[j]
                        if p1 != p2
                            println("Error: particles differ. ")
                            @show p1
                            @show p2
                            differ = true
                            break
                        end
                    end
                    differ && break
                end
                if sum(p -> p.index, ci.particles[1:ci.n_particles]) !=
                   sum(p -> p.index, cj.particles[1:cj.n_particles])
                    println("Error: particles of cells differ")
                    differ = true
                end
                if differ
                    println("Error: on cell with linear_index: $(ci.linear_index)")
                end
            end
        end
    end
    println(" Checked $ncheck of $(cl1.n_cells_with_particles) cells.")
    println(" Checked $nparticle_check of $(cl1.n_particles) particles.")
end

function check_cl(cl, box; mark)
    for i in 1:cl.n_cells_with_particles-1
        cellᵢ = cl.cells[i]
        for j in i+1:cl.n_cells_with_particles
            cellⱼ = cl.cells[j]
            for ip in 1:cellᵢ.n_particles
                xᵢ = cellᵢ.particles[ip]
                for jp in 1:cellⱼ.n_particles
                    xⱼ = cellⱼ.particles[jp]
                    if norm(xᵢ.coordinates - xⱼ.coordinates) ≈ 0
                        println("duplicate found:")
                        @show mark
                        @show Threads.threadid()
                        @show i, j, cellᵢ.linear_index, cellⱼ.linear_index
                        @show cellᵢ.n_particles
                        @show cellⱼ.n_particles
                        @show xᵢ
                        @show xⱼ
                        @show particle_cell(xᵢ.coordinates, box)
                        @show particle_cell(xⱼ.coordinates, box)
                        @show cell_linear_index(box.nc, particle_cell(xᵢ.coordinates, box))
                        @show cell_linear_index(box.nc, particle_cell(xⱼ.coordinates, box))
                        error()
                    end
                end
            end
        end
    end
end

# Create a set of random particles with an atomic density simiilar to that
# of water at room temperature and pressure.
function xatomic(n=100_000; cutoff=12.0, lcell=1, type=:vec)
    atomic_density_of_water = 0.1 # atoms/ Angstrom^3
    vol = n / atomic_density_of_water
    l = vol^(1 / 3)
    if type == :mat
        x = l * rand(3, n)
    elseif type == :vec
        x = [l * rand(SVector{3,Float64}) for _ in 1:n]
    end
    box = Box([l, l, l], cutoff, lcell=lcell)
    return x, box
end

# Create a set of  random particles with the typical density of galaxies
function xgalactic(n=8_000_000; cutoff=5.0, lcell=1, type=:vec)
    galactic_density = 12 # galaxyes per cubic mega-parsecs (https://doi.org/10.1093/mnras/106.2.121)
    vol = n / galactic_density
    l = vol^(1 / 3) # mega-parsecs
    if type == :mat
        x = l * rand(3, n)
    elseif type == :vec
        x = [l * rand(SVector{3,Float64}) for _ in 1:n]
    end
    box = Box([l, l, l], cutoff, lcell=lcell)
    return x, box
end



function test_pathological(; nmax=1000, npoints=2)

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

    n = 0
    while n < nmax
        x = rand(SVector{2,Float64}, npoints)
        matrix = rand(matrices)
        lcell = rand(1:5)
        cutoff = 0.2
        box = Box(matrix, cutoff; lcell=lcell)
        naive = CellListMap.map_naive!((x, y, i, j, d2, out) -> g(i, j, d2, box.cutoff, out), 0, x, box)
        cl = CellList(x, box)
        clmap = map_pairwise((x, y, i, j, d2, out) -> g(i, j, d2, box.cutoff, out), 0, box, cl)
        if naive != clmap
            return x, box
        end
        n += 1
    end
    return nothing, nothing
end

@testitem "test_pathological2D" begin
    @test CellListMap.test_pathological() == (nothing, nothing)
end



