"""

```julia
map_naive!(f,output,x,box)
```

Internal function or structure - interface may change.

# Extended help

Function that uses the naive pairwise mapping algorithm, for testing.

"""
function map_naive!(f, output, x, box::Box)
    @unpack unit_cell, cutoff_sq = box
    for i in 1:length(x) - 1
        xᵢ = x[i]
        for j in i + 1:length(x)
            xⱼ = wrap_relative_to(x[j], xᵢ, box)
            d2 = norm_sqr(xᵢ - xⱼ)
            if d2 <= cutoff_sq
                output = f(xᵢ, xⱼ, i, j, d2, output)
            end
        end
    end
    return output
end

function map_naive!(f, output, x, y, box::Box)
    @unpack unit_cell, cutoff_sq = box
    for i in 1:length(x)
        xᵢ = x[i]
        for j in 1:length(y)
            yⱼ = wrap_relative_to(y[j], xᵢ, box)
            d2 = norm_sqr(xᵢ - yⱼ)
            if d2 <= cutoff_sq
                output = f(xᵢ, yⱼ, i, j, d2, output)
            end
        end
    end
    return output
end

"""

```
view_celllist_particles(cl::CellList)
```

Internal function or structure - interface may change.

# Extended help

Auxiliary function to view the particles of a computing box, including images created
for computing purposes.

## Example

```julia-repl
julia> box = Box([ 100 50; 50 100 ],10);

julia> x = [ box.unit_cell_max .* rand(SVector{2,Float64}) for i in 1:1000 ];

julia> cl = CellList(x,box);

julia> p = CellListMap.view_celllist_particles(cl);

julia> using Plots

julia> scatter(Tuple.(p),label=nothing,xlims=(-10,180),ylims=(-10,180))

```

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

test_map2(box,cl;parallel=true) = 
    map_pairwise!(
        (x, y, i, j, d2, s) -> begin
            s += 1
        end,
        0, box, cl, parallel=parallel
    )

function missing_pair(x1,x2)
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

test_map(box,cl;parallel=true) = map_pairwise!((x, y, i, j, d2, s) -> s += 1/d2, 0., box, cl, parallel=parallel)
test_naive(x,box) = CellListMap.map_naive!((x, y, i, j, d2, s) -> s += 1/d2, 0., x, box)
function simple_test(x,box;parallel=true)
    cl = CellList(x,box,parallel=parallel)
    r_naive = test_naive(x,box)
    r_map = test_map(box,cl,parallel=parallel)
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
                @set! unit_cell_matrix[i,i] = 1 + 10*rand()
            end
        else
            for i in 1:N, j in 1:N
                @set! unit_cell_matrix[i,j] = 1 + 10*rand()
            end
        end
        cutoff = 1 + rand()
        try 
            if !check_unit_cell(unit_cell_matrix, cutoff, printerr=false)
                continue
            end
        catch
            @show unit_cell_matrix, cutoff
            return false
        end
        if UnitCellType == OrthorhombicCell
            box = Box([unit_cell_matrix[i,i] for i in 1:N], cutoff, lcell=lcell)    
        else
            box = Box(unit_cell_matrix, cutoff, lcell=lcell)    
        end
        if prod(box.nc) > 100000
            continue
        end
        show_progress && (next!(p))
        x = 10 * rand(SVector{N,Float64}, M)
        for i in eachindex(x)
            x[i] = x[i] .- 50  
        end
        cl = CellList(x, box, parallel=parallel)
        test = test_map(box, cl, parallel=parallel) 
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
    return true, x, box
end

function drawbox(box::Box{UnitCellType,2}) where UnitCellType
    S = SVector{2,Float64}
    m = box.unit_cell.matrix
    x = [
        S(0., 0.),
        S(m[:,1]),
        S(m[:,1] + m[:,2]),
        S(m[:,2]),
        S(0., 0.)
    ]
    return x
end

function drawbox(box::Box{UnitCellType,3}) where UnitCellType
    S = SVector{3,Float64}
    m = box.unit_cell.matrix
    x = [
        S(0., 0., 0.),
        S(m[:,1]),
        S(m[:,1] + m[:,2]),
        S(m[:,2]),
        S(0., 0., 0.),
        S(m[:,1]),
        S(m[:,1] + m[:,3]),
        S(m[:,3]),
        S(0., 0., 0.),
        S(m[:,2]),
        S(m[:,2] + m[:,3]),
        S(m[:,2]),
        S(0., 0., 0.),
        S(m[:,1]),
        S(m[:,1] + m[:,2]),
        S(m[:,1] + m[:,2]) + m[:,3],
        S(m[:,1] + m[:,3]),
        S(m[:,3]),
        S(m[:,3]) + m[:,2],
        S(m[:,1] + m[:,2]) + m[:,3]
    ]
    return x
end

"""

```
draw_computing_cell(x,box::Box{UnitCellType,2}) where UnitCellType
```

```
draw_computing_cell(cl::CellList,box::Box{UnitCellType,2},x) where UnitCellType
```

This function creates a plot of the computing cell, in two dimensions.

"""
function draw_computing_cell(x, box::Box{UnitCellType,2};parallel=true) where UnitCellType
    cl = CellList(x, box, parallel=parallel)
    return draw_computing_cell(cl,box,x) 
end
function draw_computing_cell(cl::CellList, box::Box{UnitCellType,2}, x) where UnitCellType
    box_points = drawbox(box)
    p = view_celllist_particles(cl)
    plt = Main.plot()
    Main.plot!(plt, Tuple.(box_points), label=:none)
    Main.scatter!(plt, Tuple.(p), label=:none, markeralpha=0.3)
    Main.scatter!(plt, Tuple.(wrap_to_first.(x, Ref(box))), label=:none)
    xmin = minimum(el[1] for el in p) - 3 * box.cell_size[1]
    xmax = maximum(el[1] for el in p) + 3 * box.cell_size[1]
    ymin = minimum(el[2] for el in p) - 3 * box.cell_size[2]
    ymax = maximum(el[2] for el in p) + 3 * box.cell_size[2]
    Main.plot!(plt,
        aspect_ratio=1,framestyle=:box,xrotation=60,
        xlims=(xmin, xmax),
        ylims=(ymin, ymax),
        xticks=(round.(digits=3, xmin:box.cell_size[1]:xmax)),
        yticks=(round.(digits=3, ymin:box.cell_size[2]:ymax)),
    )
    return plt
end

"""

```
draw_computing_cell(x,box::Box{UnitCellType,3}) where UnitCellType
```

This function creates a plot of the computing cell, in three dimensions.

"""
function draw_computing_cell(x, box::Box{UnitCellType,3}; parallel=true) where UnitCellType
    cl = CellList(x, box, parallel=parallel)
    box_points = drawbox(box)
    p = view_celllist_particles(cl)
    plt = Main.plot()
    Main.plot!(plt, Tuple.(box_points), label=:none)
    Main.scatter!(plt, Tuple.(p), label=:none, markeralpha=0.3)
    Main.scatter!(plt, Tuple.(wrap_to_first.(x, Ref(box))), label=:none, markeralpha=0.3)
    lims = Vector{Float64}[] 
    for i in 1:3
        push!(lims, [ -2 * box.cell_size[i], box.nc[i] + 2 * box.cell_size[i] ])
    end
    Main.plot!(plt,
        aspect_ratio=1,framestyle=:box,xrotation=60,yrotation=-70,zrotation=0,
        xlims=lims[1],
        ylims=lims[2],
        zlims=lims[3],
        xticks=(round.(digits=3, lims[1][1]:box.cell_size[1]:lims[1][2])),
        yticks=(round.(digits=3, lims[2][1]:box.cell_size[2]:lims[2][2])),
        zticks=(round.(digits=3, lims[3][1]:box.cell_size[3]:lims[3][2])),
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

function check_cl(cl,box;mark)
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
                        @show particle_cell(xᵢ.coordinates,box)
                        @show particle_cell(xⱼ.coordinates,box)
                        @show cell_linear_index(box.nc,particle_cell(xᵢ.coordinates,box))
                        @show cell_linear_index(box.nc,particle_cell(xⱼ.coordinates,box))
                        error()
                    end
                end
            end
        end
    end
end

# Create a set of random particles with an atomic density simiilar to that
# of water at room temperature and pressure.
function xatomic(n=100_000;cutoff=12.,lcell=1)
    atomic_density_of_water = 0.1 # atoms/ Angstrom^3
    vol = n/atomic_density_of_water
    l = vol^(1/3)
    x = l*rand(3,n)
    box = Box([l,l,l],cutoff,lcell=lcell)
    return x, box
end

# Create a set of  random particles with the typical density of galaxies
function xgalactic(n=8_000_000;cutoff=5.,lcell=1)
    galactic_density = 12 # galaxyes per cubic mega-parsecs (https://doi.org/10.1093/mnras/106.2.121)
    vol = n/galactic_density
    l = vol^(1/3) # mega-parsecs
    x = l*rand(3,n)
    box = Box([l,l,l],cutoff,lcell=lcell)
    return x, box
end






