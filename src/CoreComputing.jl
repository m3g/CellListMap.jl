#
# Parallel thread spliiter
#
splitter(first,n) = first:nthreads():n

"""

```
reduce(output::Number, output_threaded::Vector{<:Number})
```

Functions to reduce the output of common options (vectors of numbers 
and vectors of vectors). This function can be overloaded by custom
reduction methods. It always must both receive the `output` variable
as a parameter, and return it at the end.

"""
reduce(output::Number, output_threaded::Vector{<:Number}) = sum(output_threaded)
function reduce(output::AbstractVector, output_threaded::AbstractVector{<:AbstractVector}) 
  for i in 1:nthreads()
    @. output += output_threaded[i]
  end
  return output
end

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, cl::CellList{N,T}; 
  show_progress::Bool=false
) where {F,N,T}
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
  f::F1, output, box::Box, cl::CellList{N,T};
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress::Bool=false
) where {F1,F2,N,T}
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

"""

```
partition!(x::AbstractVector,by)
```

Function that reorders `x` vector by putting in the first positions the
elements with values satisfying `by(el)`. Returns the number of elements
that satisfy the condition.

"""
function partition!(x::AbstractVector,by)
  iswap = 1
  @inbounds for i in eachindex(x)
    if by(x[i])
      if iswap != i
        x[iswap], x[i] = x[i], x[iswap]
      end
      iswap += 1
    end
  end
  return iswap - 1
end


"""

```
project_particles(jc_cartesian,cl,Δc,icell,box)
```

Projects all particles of the cell of cartesian coordinates
`jc_cartesian` into the unnitary vector `Δc` with direction
connecting `icell.center` and the center of cell `jc_cartesian`. 
Returns a view of the array of projected particles, with length
equal to the number of particles of cell `jc_cartesian`. 

"""
function project_particles(jc_cartesian,cl::CellList,Δc,icell,box)
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
    xreal = pⱼ.real
    projected_particles[jp] = ProjectedParticle(j_orig,xproj,xpⱼ,xreal) 
    pⱼ = cl.np[j]
    j = pⱼ.index
  end
  pp = @view(projected_particles[1:npcell])
  return pp
end

function inner_loop!(
  f,box::Box{TriclinicCell},icell,
  cl::CellList{N,T},
  output
) where {N,T}
  @unpack cutoff, cutoff_sq, nc = box
  cell = cl.cwp[icell]

  for neighbour_cell in neighbour_cells_all(box)
    jc_cartesian = cell.cartesian + neighbour_cell
    jc = cell_linear_index(nc,jc_cartesian)

    # Vector connecting cell centers
    if jc == cell.icell
      Δc = SVector{N,T}(ntuple(i->1,N))
    else
      Δc = cell_center(jc_cartesian,box) - cell.center 
    end
    Δc = Δc/norm(Δc)
    pp = project_particles(jc_cartesian,cl,Δc,cell,box)

    pᵢ = cl.fp[cell.icell]
    i = pᵢ.index
    while i > 0
      if !pᵢ.real
        pᵢ = cl.np[pᵢ.index]
        i = pᵢ.index
        continue
      end
      xpᵢ = pᵢ.coordinates
      i_orig = pᵢ.index_original

      # Partition pp array according to the current projections. This
      # is faster than it looks, because after the first partitioning, the
      # array will be almost correct, and essentially the algorithm
      # will only run over most elements. Avoiding this partitioning or
      # trying to sort the array before always resulted to be much more
      # expensive. 
      xproj = dot(xpᵢ-cell.center,Δc)
      n = partition!(pp, el -> el.xproj - xproj <= cutoff)

      for j in 1:n
        pⱼ = pp[j]
        j_orig = pⱼ.index_original
        if i_orig < j_orig
          xpⱼ = pⱼ.coordinates
          d2 = norm_sqr(xpᵢ - xpⱼ)
          if d2 <= cutoff_sq
            output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
          end
        end
      end
      pᵢ = cl.np[pᵢ.index]
      i = pᵢ.index
    end
  end

  return output
end

#
# Inner loop for Orthorhombic cells is faster because we can guarantee that
# there are not repeated computations even if running over half of the cells.
#
function inner_loop!(
  f,box::Box{OrthorhombicCell},icell,
  cl::CellList{N,T},
  output
) where {N,T}
  @unpack cutoff_sq = box
  cell = cl.cwp[icell]

  # loop over list of non-repeated particles of cell ic
  pᵢ = cl.fp[cell.icell]
  i = pᵢ.index
  while i > 0
    xpᵢ = pᵢ.coordinates
    pⱼ = cl.np[i] 
    j = pⱼ.index # here j starts from i+1
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

function cell_output!(
  f,
  box::Box{OrthorhombicCell},
  icell,
  cl::CellList{N,T},
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
      pⱼ = pp[j]
      xpⱼ = pⱼ.coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      if d2 <= cutoff_sq
        i_orig = pᵢ.index_original
        j_orig = pⱼ.index_original
        output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      end
    end
    pᵢ = cl.np[i]
    i = pᵢ.index
  end

  return output
end

#
# Serial version for cross-interaction computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, 
  cl::CellListPair{N,T}; 
  show_progress=show_progress
) where {F,N,T}
  show_progress && (p = Progress(length(cl.small),dt=0))
  for i in eachindex(cl.small)
    output = inner_loop!(f,output,i,box,cl)
    show_progress && next!(p)
  end
  return output
end

#
# Parallel version for cross-interaction computations
#
function map_pairwise_parallel!(
  f::F1, output, box::Box, 
  cl::CellListPair{N,T};
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress=show_progress
) where {F1,F2,N,T}
  show_progress && (p = Progress(length(cl.small),dt=1))
  @threads for it in 1:nthreads()
    for i in splitter(it,length(cl.small))
      output_threaded[it] = inner_loop!(f,output_threaded[it],i,box,cl) 
      show_progress && next!(p)
    end
  end 
  output = reduce(output,output_threaded)
  return output
end

#
# Inner loop of cross-interaction computations. If the two sets
# are large, this might(?) be improved by computing two sepparate
# cell lists and using projection and partitioning.
#
function inner_loop!(
  f,output,i,box,
  cl::CellListPair{N,T}
) where {N,T}
  @unpack nc, cutoff_sq = box
  xpᵢ = wrap_to_first(cl.small[i],box)
  ic = particle_cell(xpᵢ,box)
  for neighbour_cell in neighbour_cells_all(box)
    jc = cell_linear_index(nc,neighbour_cell+ic)
    pⱼ = cl.large.fp[jc]
    j = pⱼ.index
    # loop over particles of cell jc
    while j > 0
      xpⱼ = pⱼ.coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      if d2 <= cutoff_sq
        j_orig = pⱼ.index_original 
        output = f(xpᵢ,xpⱼ,i,j_orig,d2,output)
      end
      pⱼ = cl.large.np[j]
      j = pⱼ.index
    end                                   
  end
  return output
end

#
# Serial version for self-pairwise computations
#
function map_pairwise_serial!(
  f::F, output, box::Box, cl::CellListSIMD{N,T}; 
  show_progress::Bool=false
) where {F,N,T}
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
  f::F1, output, box::Box, cl::CellListSIMD{N,T};
  output_threaded=output_threaded,
  reduce::F2=reduce,
  show_progress::Bool=false
) where {F1,F2,N,T}
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
  f,box::Box{TriclinicCell},icell,
  cl::CellListSIMD{N,T},
  output
) where {N,T}
  @unpack cutoff, cutoff_sq, nc = box
  cell = cl.cwp[icell]

  for neighbour_cell in neighbour_cells_all(box)
    jc_cartesian = cell.cartesian + neighbour_cell
    jc = cell_linear_index(nc,jc_cartesian)

    # Vector connecting cell centers
    if jc == cell.icell
      Δc = SVector{N,T}(ntuple(i->1,N))
    else
      Δc = cell_center(jc_cartesian,box) - cell.center 
    end
    Δc = Δc/norm(Δc)
    pp = project_particles(jc_cartesian,cl,Δc,cell,box)

    for i in 1:cl.npcell[cell.icell]
      pᵢ = cl.list[cell.icell][i]
      (!pᵢ.real) && continue
      xpᵢ = pᵢ.coordinates
      i_orig = pᵢ.index_original

      # Partition pp array according to the current projections. This
      # is faster than it looks, because after the first partitioning, the
      # array will be almost correct, and essentially the algorithm
      # will only run over most elements. Avoiding this partitioning or
      # trying to sort the array before always resulted to be much more
      # expensive. 
      xproj = dot(xpᵢ-cell.center,Δc)
      n = partition!(pp, el -> el.xproj - xproj <= cutoff)

      for j in 1:n
        pⱼ = pp[j]
        j_orig = pⱼ.index_original
        if i_orig < j_orig
          xpⱼ = pⱼ.coordinates
          d2 = norm_sqr(xpᵢ - xpⱼ)
          if d2 <= cutoff_sq
            output = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
          end
        end
      end
      pᵢ = cl.np[pᵢ.index]
      i = pᵢ.index
    end
  end

  return output
end

#
# Inner loop for Orthorhombic cells is faster because we can guarantee that
# there are not repeated computations even if running over half of the cells.
#
function inner_loop!(
  f,box::Box{OrthorhombicCell},icell,
  cl::CellListSIMD{N,T},
  output
) where {N,T}
  @unpack cutoff_sq = box
  cell = cl.cwp[icell]

  # loop over list of non-repeated particles of cell ic
  for i in 1:cl.npcell[cell.icell]-1
    pᵢ = cl.list[cell.icell][i]
    xpᵢ = pᵢ.coordinates
    for j in i+1:cl.npcell[cell.icell]
      pⱼ = cl.list[cell.icell][j]
      xpⱼ = pⱼ.coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      within_cutoff = (d2 <= cutoff_sq)
      i_orig = pᵢ.index_original
      j_orig = pⱼ.index_original
      output_temp = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      output = (1-within_cutoff)*output + within_cutoff*output_temp
    end
  end

  for jcell in neighbour_cells(box)
    output = cell_output!(f,box,cell,cl,output,cell.cartesian+jcell)
  end

  return output
end

function cell_output!(
  f,
  box::Box{OrthorhombicCell},
  icell,
  cl::CellListSIMD{N,T},
  output,
  jc_cartesian
) where {N,T}
  @unpack cutoff, cutoff_sq = box

  # Vector connecting cell centers
  Δc = cell_center(jc_cartesian,box) - icell.center 
  Δc = Δc/norm(Δc)
  pp = project_particles(jc_cartesian,cl,Δc,icell,box)

  # Loop over particles of cell icell
  for i in 1:cl.npcell[icell.icell]
    pᵢ = cl.list[icell.icell][i]
    xpᵢ = pᵢ.coordinates
    xproj = dot(xpᵢ-icell.center,Δc)
    
    # Partition pp array according to the current projections
    n = partition!(pp, el -> el.xproj - xproj <= cutoff)

    # Compute the interactions 
    for j in 1:n
      pⱼ = pp[j]
      xpⱼ = pⱼ.coordinates
      d2 = norm_sqr(xpᵢ - xpⱼ)
      within_cutoff = (d2 <= cutoff_sq)
      i_orig = pᵢ.index_original
      j_orig = pⱼ.index_original
      output_temp = f(xpᵢ,xpⱼ,i_orig,j_orig,d2,output)
      output = (1-within_cutoff)*output + within_cutoff*output_temp
    end
  end

  return output
end

"""

```
project_particles(jc_cartesian,cl,Δc,icell,box)
```

Projects all particles of the cell of cartesian coordinates
`jc_cartesian` into the unnitary vector `Δc` with direction
connecting `icell.center` and the center of cell `jc_cartesian`. 
Returns a view of the array of projected particles, with length
equal to the number of particles of cell `jc_cartesian`. 

"""
function project_particles(jc_cartesian,cl::CellListSIMD,Δc,icell,box)
  @unpack nc = box
  @unpack npcell, list = cl
  projected_particles = cl.projected_particles[threadid()]
  jc = cell_linear_index(nc,jc_cartesian)
  for j in 1:npcell[jc]
    pⱼ = list[jc][j]
    j_orig = pⱼ.index_original
    xpⱼ = pⱼ.coordinates
    xproj = dot(xpⱼ-icell.center,Δc)
    xreal = pⱼ.real
    projected_particles[j] = ProjectedParticle(j_orig,xproj,xpⱼ,xreal) 
  end
  pp = @view(projected_particles[1:npcell[jc]])
  return pp
end