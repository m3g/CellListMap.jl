#
# Functions to deal with large and dense system, avoiding the maximum number
# of unnecessary loop iterations over non-interacting particles
#

const UnionLargeDense = Union{LargeDenseSystem,HugeDenseSystem}

"""

```
UpdateCellList!(
    x::AbstractVector{SVector{N,T}},
    box::Box,cl:CellList{LargeDenseSystem,N,T},
    parallel=true
) where {N,T}
```

Function that will update a previously allocated `CellList` structure, given new updated particle 
positions of large and dense systems.

## Example

```julia-repl
julia> box = Box([250,250,250],10);

julia> x = [ 250*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> cl = CellList(x,box);

julia> box = Box([260,260,260],10);

julia> x = [ 260*rand(SVector{3,Float64}) for i in 1:1000 ];

julia> cl = UpdateCellList!(x,box,cl); # update lists

```

"""
function UpdateCellList!(
  x::AbstractVector{SVector{N,T}},
  box::Box,
  cl::CellList{SystemType,N,T};
  parallel::Bool=true
) where {SystemType<:UnionLargeDense,N,T}
  @unpack contains_real, cwp, fp, np, npcell, projected_particles = cl

  number_of_cells = prod(box.nc)
  if number_of_cells > length(cwp) 
    number_of_cells = ceil(Int,1.2*number_of_cells) # some margin in case of box size variations
    resize!(contains_real,number_of_cells)
    resize!(cwp,number_of_cells)
    resize!(fp,number_of_cells)
    resize!(npcell,number_of_cells)
  end

  cl.ncwp[1] = 0
  if parallel
    @threads for i in eachindex(cwp)
      contains_real[i] = false
      cwp[i] = zero(Cell{N,T})
      fp[i] = zero(AtomWithIndex{N,T})
      npcell[i] = 0
    end
    @threads for i in eachindex(np)
      np[i] = zero(AtomWithIndex{N,T})
    end
  else
    fill!(contains_real,false)
    fill!(cwp,zero(Cell{N,T}))
    fill!(fp,zero(AtomWithIndex{N,T}))
    fill!(np,zero(AtomWithIndex{N,T}))
    fill!(npcell,0)
  end

  #
  # The following part cannot be *easily* paralelized, because 
  # there is concurrency on the construction of the cell lists
  #

  #
  # Add virtual particles to edge cells
  #
  for (ip,particle) in pairs(x)
    p = wrap_to_first(particle,box.unit_cell)
    cl = replicate_particle!(ip,p,box,cl)
  end
  #
  # Add true particles, such that the first particle of each cell is
  # always a true particle
  #
  for (ip,particle) in pairs(x)
    p = wrap_to_first(particle,box.unit_cell)
    cl = add_particle_to_celllist!(ip,p,box,cl) 
  end

  maximum_npcell = maximum(npcell)
  if maximum_npcell > length(projected_particles[1])
    for i in 1:nthreads()
      resize!(projected_particles[i],ceil(Int,1.2*maximum_npcell))
    end
  end

  return cl
end

"""

Set one index of a cell list

"""
function add_particle_to_celllist!(
  ip,
  x::SVector{N,T},
  box,
  cl::CellList{SystemType,N,T};
  real_particle::Bool=true
) where {SystemType<:UnionLargeDense,N,T}
  @unpack contains_real, ncp, ncwp, cwp, fp, np, npcell = cl
  ncp[1] += 1
  icell_cartesian = particle_cell(x,box)
  icell = cell_linear_index(box.nc,icell_cartesian)
  #
  # Cells starting with real particles are annotated to be run over
  #
  if real_particle && (!contains_real[icell])
    contains_real[icell] = true
    ncwp[1] += 1
    cwp[ncwp[1]] = Cell{N,T}(
      icell,
      icell_cartesian,
      cell_center(icell_cartesian,box)
    )
  end
  if fp[icell].index == 0
    npcell[icell] = 1
  else
    npcell[icell] += 1
  end
  if ncp[1] > length(np) 
    old_length = length(np)
    resize!(np,ceil(Int,1.2*old_length))
    for i in old_length+1:length(np)
      np[i] = zero(AtomWithIndex{N,T}) 
    end
  end
  np[ncp[1]] = fp[icell]
  fp[icell] = AtomWithIndex(ncp[1],ip,x) 
  return cl
end

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, cl::CellList{SystemType,N,T}; 
  show_progress::Bool=false
) where {F,SystemType<:UnionLargeDense,N,T}
  show_progress && (p = Progress(cl.ncwp[1],dt=1))
  for icell in 1:cl.ncwp[1]
    output = inner_loop!(f,box,icell,cl,output) 
    show_progress && next!(p)
  end
  return output
end

#
# Parallel version for self-pairwise computations
#
function map_pairwise_parallel!(
  f::F1, output, box::Box, cl::CellList{SystemType,N,T};
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress::Bool=false
) where {F1,F2,SystemType<:UnionLargeDense,N,T}
  show_progress && (p = Progress(cl.ncwp[1],dt=1))
  @threads for it in 1:nthreads() 
    for icell in splitter(it,cl.ncwp[1])
      output_threaded[it] = inner_loop!(f,box,icell,cl,output_threaded[it]) 
      show_progress && next!(p)
    end
  end 
  output = reduce(output,output_threaded)
  return output
end

function inner_loop!(
  f,box,icell,
  cl::CellList{SystemType,N,T},
  output
) where {SystemType<:UnionLargeDense,N,T}
  @unpack cutoff_sq = box
  cell = cl.cwp[icell]

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[cell.icell]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    pⱼ = cl.np[i] 
    j = pⱼ.index
    while j > 0
      xpⱼ = pⱼ.coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      if d2 <= cutoff_sq
        i_orig = pᵢ.index_original
        j_orig = pⱼ.index_original
          output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      end
      pⱼ = cl.np[pⱼ.index]
      j = pⱼ.index
    end
    pᵢ = cl.np[pᵢ.index]
    i = pᵢ.index
  end

  for jcell in neighbour_cells(box)
    output = cell_output!(f,box,cell,cl,output,cell.cartesian+jcell)
  end

  return output
end

function project_particles(jc_cartesian,cl,Δc,icell,box)
  @unpack nc = box
  @unpack fp = cl
  projected_particles = cl.projected_particles[threadid()]
  jc = cell_linear_index(nc,jc_cartesian)
  pⱼ = fp[jc]
  npcell = cl.npcell[jc]
  j = pⱼ.index
  for jp in 1:npcell
    j_orig = pⱼ.index_original
    xpⱼ = pⱼ.coordinates
    xproj = dot(xpⱼ-icell.center,Δc)
    projected_particles[jp] = ProjectedParticle(j_orig,xproj,xpⱼ) 
    pⱼ = cl.np[j]
    j = pⱼ.index
  end
  pp = @view(projected_particles[1:npcell])
  return pp
end

#
# loops over the particles of a neighbour cell, for HugeDenseSystem
#
# all efforts are made to not run anything on unnecessary pairs 
# of particles, this includes sorting the projections
#
function cell_output!(
  f,
  box,
  icell,
  cl::CellList{HugeDenseSystem,N,T},
  output,
  jc_cartesian
) where {N,T}
  @unpack cutoff, cutoff_sq = box

  # Project particles into vector connection cell centers
  Δc = cell_center(jc_cartesian,box) - icell.center 
  Δc = Δc/norm(Δc)
  pp = project_particles(jc_cartesian,cl,Δc,icell,box)

  # Sort particles according to projection norm
  sort!(pp, by=el->el.xproj,alg=InsertionSort)

  # Loop over particles of cell icell
  pᵢ = cl.fp[icell.icell]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    xproj = dot(xpᵢ-icell.center,Δc)
    for j in eachindex(pp)
      xproj - pp[j].xproj > cutoff && break
      xpⱼ = pp[j].coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      if d2 <= cutoff_sq
        i_orig = pᵢ.index_original
        j_orig = pp[j].index_original
        output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      end
    end
    pᵢ = cl.np[i]
    i = pᵢ.index
  end

  return output
end

#
# loops over the particles of a neighbour cell, for LargeDenseSystem
# not large enough such that sorting the projections is useful
#
function cell_output!(
  f,
  box,
  icell,
  cl::CellList{LargeDenseSystem,N,T},
  output,
  jc_cartesian
) where {N,T}
  @unpack cutoff, cutoff_sq = box

  # Vector connecting cell centers
  Δc = cell_center(jc_cartesian,box) - icell.center 
  Δc = Δc/norm(Δc)
  pp = project_particles(jc_cartesian,cl,Δc,icell,box)

  # Loop over particles of cell icell
  pᵢ = cl.fp[icell.icell]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    xproj = dot(xpᵢ-icell.center,Δc)

    # Partition pp array according to the current projections
    n = partition!(pp, el -> el.xproj - xproj <= cutoff)

    # Compute the interactions 
    for j in 1:n
      xpⱼ = pp[j].coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      if d2 <= cutoff_sq
        i_orig = pᵢ.index_original
        j_orig = pp[j].index_original
        output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      end
    end
    pᵢ = cl.np[i]
    i = pᵢ.index
  end

  return output
end
