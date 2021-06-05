module CellListMap

using Base.Threads
using Parameters
using StaticArrays
using DocStringExtensions

export LinkedLists, Box, initcells!
export map_pairwise

"""

```
LinkedLists(N)
```

Structure that contains the vectors storing the first particle and next particle of the linked
cells. To be initialized with the number of particles.

## Example

```julia-repl
julia> lc = LinkedLists(100_000)

julia> lc = LinkedLists(100_000)
LinkedLists{100000}
  firstatom: Array{Int64}((100000,)) [0, 0, 0, 0, 0, 0, 0, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  nextatom: Array{Int64}((100000,)) [0, 0, 0, 0, 0, 0, 0, 0, 0, 0  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

```

"""  
@with_kw struct LinkedLists{N}
  firstatom::Vector{Int}
  nextatom::Vector{Int}
end
function LinkedLists(N::Int) 
  return LinkedLists{N}(
    zeros(Int,N), # firstatom (actual required size is nc1*nc2*nc*3
    zeros(Int,N), # nextatom
  )
end

"""

$(TYPEDEF)

Structure that contains some data required to compute the linked cells. To
be initialized with the box size and cutoff. An optional parameter `lcell` 
can be provided as the last argument to define the cell size relative to the
cutoff (default 2, meaning half of the cutoff).

## Example

```julia-repl
julia> sides = [250,250,250];

julia> cutoff = 10.;

julia> box = Box(sides,cutoff)

julia> box = Box(sides,cutoff)
Box{3, Float64}([250.0, 250.0, 250.0], [50, 50, 50], [5.0, 5.0, 5.0], 2, 10.0)

```


"""
@with_kw struct Box{N,T<:Float64}
  sides::SVector{N,T}
  nc::SVector{N,Int}
  l::SVector{N,T}
  lcell::Int
  cutoff::T
end
function Box(sides::SVector{N,T}, cutoff, lcell::Int=2) where {N,T<:Float64}
  # Compute the number of cells in each dimension
  nc = SVector{N,Int}(max.(1,trunc.(Int,sides/(cutoff/lcell))))
  l = SVector{N,T}(sides ./ nc)
  return Box{N,T}(sides,nc,l,lcell,cutoff)
end
function Box(sides::AbstractVector{T}, cutoff, lcell::Int=2) where T 
  N =length(sides)
  Box(SVector{N,Float64}(sides), cutoff, lcell)
end

"""

```
sq_distance(x,y)
```

Function to compute squared Euclidean distances between two n-dimensional vectors.

"""
@inline function sq_distance(x::AbstractVector{T}, y::AbstractVector{T}) where T
  @assert length(x) == length(y)
  d = zero(T)
  @inbounds for i in eachindex(x)
    d += (x[i]-y[i])^2
  end
  return d
end

"""

```
distance(x,y)
```

Function to compute Euclidean distances between two n-dimensional vectors.

"""
@inline distance(x::AbstractVector{T}, y::AbstractVector{T}) where T =
  sqrt(sq_distance(x,y))

"""

```
particle_cell(x::AbstractVector{T}, box::Box{N,T}) where {N,T}
```

Returns the coordinates of the cell to which a particle belongs, given its coordinates
and the sides of the periodic box (for arbitrary dimension N).

"""
function particle_cell(x::AbstractVector{T}, box::Box{N,T}) where {N,T}
  # Wrap to origin
  xwrapped = wrapone(x,box.sides)
  return ntuple(i -> floor(Int,(xwrapped[i]+box.sides[i]/2)/box.l[i])+1, N)
end
"""

```
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N}
```

Given the linear index of the cell in the cell list, returns the cartesian indexes
of the cell (for arbitrary dimension N).


"""
cell_cartesian_indices(nc::SVector{N,Int}, i1D) where {N} = 
  CartesianIndices(ntuple(i -> nc[i],N))[i1D]

"""
```
icell1D(nc::SVector{N,Int}, i, j, k) where N
```
Returns the index of the cell, in the 1D representation, from its cartesian coordinates.

"""
cell_linear_index(nc::SVector{N,Int}, indexes...) where {N} =
  LinearIndices(ntuple(i -> nc[i],N))[ntuple(i->indexes[i],N)...]

"""

```
function wrap!(x::AbstractVector{T}, sides::T, center::T) where T <: AbstractVector
```

Functions that wrap the coordinates They modify the coordinates of the input vector.  
Wrap to a given center of coordinates

"""
function wrap!(x::AbstractVector{T}, 
               sides::T, 
               center::T) where T <: AbstractVector
  for i in eachindex(x)
    x[i] = wrapone(x[i],sides,center)
  end
  return nothing
end

@inline function wrapone(x::T, sides::T, center::T) where T <: AbstractVector
  s = @. (x-center)%sides
  s = @. wrapx(s,sides) + center
  return s
end

@inline function wrapx(x,s)
  if x > s/2
    x = x - s
  elseif x < -s/2
    x = x + s
  end
  return x
end

"""

```
wrap!(x::AbstractVector{T}, sides::T) where T <: AbstractVector
```
Wrap to origin (slightly cheaper).

"""
function wrap!(x::AbstractVector{T}, sides::T) where T <: AbstractVector
  for i in eachindex(x)
    x[i] = wrapone(x[i],sides)
  end
  return nothing
end

@inline function wrapone(x::T, sides::T) where T <: AbstractVector
  s = @. x%sides
  s = @. wrapx(s,sides)
  return s
end

"""

```
wrap_cell(nc::SVector{N,T}, indexes...) where {N,T}
```

Given the `N` indexes of a cell, return the periodic cell which correspondst to
it, if the cell is outside the main box.

"""
function wrap_cell(nc::SVector{N,T}, indexes::Int...) where {N,T}
  cell_indexes = ntuple(N) do i
    ind = indexes[i]
    if ind < 1
      ind = nc[i] + ind
    elseif ind > nc[i]
      ind = ind - nc[i]
    end
    return ind
  end
  return cell_indexes
end

"""

```
initcells!(x::AbstractVector{T}, box::Box, lc::LinkedLists) where T
```
Function that initializes the linked cells by computing to each cell each atom
belongs and filling up the firstatom and nexatom arrays.
Modifies the data in the lc structure

## Example

```julia-repl
julia> lc = LinkedLists(100_000);

julia> box = Box([250,250,250],10);

julia> x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:100_000 ];

julia> initcells!(x,box,lc)

julia> lc.firstatom
125000-element Vector{Int64}:
 51137
 62998
     ⋮
 41753
 53909

```

"""
function initcells!(x::AbstractVector{<:AbstractVector}, box::Box{N,T}, lc::LinkedLists) where {N,T}

  # Count the number of boxes and checks if there is a problem with dimensions
  nboxes = prod(box.nc)
  if length(lc.firstatom) < nboxes
    resize!(lc.firstatom,nboxes)
  end

  # Reset arrays
  for i in 1:nboxes
    lc.firstatom[i] = 0
  end
  @. lc.nextatom = 0
 
  # Initialize cell, firstatom and nexatom
  for iat in eachindex(x)
    ic, jc, kc = particle_cell(x[iat],box)
    icell = cell_linear_index(box.nc,ic,jc,kc)
    lc.nextatom[iat] = lc.firstatom[icell]
    lc.firstatom[icell] = iat
  end

  return nothing
end

"""

```
map_pairwise(f::Function,output,x::AbstractVector,box::Box{N,T},lc::LinkedLists) where {N,T}
```

This function will run over every pair of particles which are closer than `box.cutoff` and compute
the Euclidean distance between the particles, considering the periodic boundary conditions given
in the `Box` structure. If the distance is smaller than the cutoff, a function `f` of the coordinates
of the two particles will be computed. 

This function `f` receives six arguments as input: 
```
f(x,y,i,j,d2,output)
```
Which are the coordinates of one particle, the coordinates of the second particle, the index of the first particle, the index of the second particle, the squared distance between them, and the `output` variable. It has also to return the same `output` variable. Thus, `f` may or not mutate `output`, but in either case it must return it. With that, it is possible to compute an average property of the distance of the particles or, for example, build a histogram. The squared distance `d2` is computed internally for comparison with the `cutoff`, and is passed to the `f` because many times it is used for the desired computation. 

## Example

Computing the mean difference in `x` position between random particles, remembering the the number of pairs of `n` particles is `n(n-1)/2`. The function does not use the indices or the distance, such that we remove them from the parameters by using a closure.

```julia-repl
julia> n = 100_000;

julia> lc = LinkedLists(n);

julia> box = Box([250,250,250],10);

julia> x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:n ];

julia> initcells!(x,box,lc)

julia> f(x,y,sum_dx) = sum_dx + x[1] - y[1] 

julia> normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

julia> avg_dx = normalization * map_parwise((x,y,i,j,d2,sum_dx) -> f(x,y,sum_dx), 0.0, x, box, lc)

```

Computing the histogram of the distances between particles (considering the same particles as in the above example). Again,
the function does not use the indices, but uses the distance, which are removed from the function call using a closure:

```
julia> function build_histogram!(x,y,d2,hist)
         d = sqrt(d2)
         ibin = floor(Int,d) + 1
         hist[ibin] += 1
         return hist
       end;

julia> hist = zeros(Int,10);

julia> normalization = N / (N*(N-1)/2) # (number of particles) / (number of pairs)

julia> hist = normalization * map_pairwise((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),hist,x,box,lc)

```

In this test we compute the "gravitational potential", pretending that each particle
has a different mass. In this case, the closure is used to pass the masses to the
function that computes the potential.

```julia
# masses
mass = rand(N)

# Function to be evalulated for each pair: build distance histogram
function potential(x,y,i,j,d2,u,mass)
  d = sqrt(d2)
  u = u - 9.8*mass[i]*mass[j]/d
  return u
end

# Run pairwise computation
u = map_pairwise((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc)
```

The example above can be run with `CellLists.test3()`.


"""
map_pairwise(
  f::Function, output, 
  x::AbstractVector,
  box::Box{N,T}, lc::LinkedLists; reduce::Function=reduce,
) where {N,T} = map_pairwise(f,output,x,x,box,lc;self=true,reduce=reduce)

"""

``
map_pairwise( f::Function, output, 
  x::AbstractVector, y::AbstractVector, 
  box::Box{N,T}, lc::LinkedLists
) where {N,T}
```

The same as the function `map_pairwise` that receives a single vector `x`, but to compute interactions between two disjoint sets of particles `x` and `y`. 

"""
function map_pairwise(
  f::Function, output, 
  x::AbstractVector, y::AbstractVector, 
  box::Box{N,T}, lc::LinkedLists; 
  self=false, reduce::Function=reduce,
) where {N,T}

  @unpack sides, nc, lcell, cutoff = box
  output_threaded = [ output ]
  for i in 2:nthreads()
    push!(output_threaded,deepcopy(output))
  end
  cutoff2 = cutoff^2
  @threads for i in eachindex(x)
    it = threadid()
    xᵢ = x[i]
    # Check the cell of this atom
    ipc, jpc, kpc = particle_cell(xᵢ,box)
    # Loop over vicinal cells to compute distances to solvent atoms, and
    # add data to dc structure (includes current cell)
    @inbounds for ic in ipc-lcell:ipc+lcell
      for jc in jpc-lcell:jpc+lcell
        for kc in kpc-lcell:kpc+lcell
          # Wrap cell if needed
          iw, jw, kw = wrap_cell(nc,ic,jc,kc)
          # get linear index of this cell
          icell = cell_linear_index(nc,iw,jw,kw) 
          # cycle over the atoms of this cell
          j = lc.firstatom[icell]
          while j > 0
            # skip same particle and repeated
            if self && j >= i 
              j = lc.nextatom[j]
              continue
            end
            # Wrap particle j relative to particle xᵢ
            yⱼ = wrapone(y[j],sides,xᵢ)
            d2 = sq_distance(xᵢ,yⱼ)
            if d2 <= cutoff2
              output_threaded[it] = f(xᵢ,yⱼ,i,j,d2,output_threaded[it])
            end
            j = lc.nextatom[j]
          end
        end
      end
    end
  end
  output = reduce(output_threaded)
  return output
end

#
# Functions to reduce the output
#
function reduce(output_threaded::Vector{<:Number}) 
  output = output_threaded[1]
  @inbounds for i in 2:nthreads()
    output = output + output_threaded[i]
  end
  return output
end

function reduce(output_threaded::Vector{<:AbstractVector}) 
  output = output_threaded[1]
  @inbounds for i in 2:nthreads()
    output .= output .+ output_threaded[i]
  end
  return output
end

#
# In this test we compute the average displacement of the x coordinates of the atoms
# Expected to be nearly zero in average
#              
function test1(N=100_000)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists
  lc = LinkedLists(N)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]

  # Initializing linked cells with these positions
  initcells!(x,box,lc)  

  # Function to be evalulated for each pair: sum of displacements on x
  f(x,y,avg_dx) = avg_dx + x[1] - y[1]

  avg_dx = (N/(N*(N-1)/2)) * map_pairwise((x,y,i,j,d2,avg_dx) -> f(x,y,avg_dx),0.,x,box,lc)
  return avg_dx

end

#
# Testing
#

#
# In this test we compute the histogram of distances, expected to follow the
# function f(f) = ρ(4/3)π(r[i+1]^3 - r[i]^3) with ρ being the density of the system.
#
function test2(N=100_000)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists
  lc = LinkedLists(N)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]

  # Initializing linked cells with these positions
  initcells!(x,box,lc)  

  # Function to be evalulated for each pair: build distance histogram
  function build_histogram!(x,y,d2,hist) 
    d = sqrt(d2)
    ibin = floor(Int,d) + 1
    hist[ibin] += 1
    return hist
  end

  # Run pairwise computation
  hist = zeros(Int,10)
  hist = map_pairwise((x,y,i,j,d2,hist) -> build_histogram!(x,y,d2,hist),hist,x,box,lc)
  return (N/(N*(N-1)/2)) * hist

end

#
# In this test we compute the "gravitational potential", pretending that each particle
# has a different mass. In this case, the closure is used to pass the masses to the
# function that computes the potential.
#
function test3(N=100_000)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists
  lc = LinkedLists(N)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]

  # masses
  mass = rand(N)

  # Initializing linked cells with these positions
  initcells!(x,box,lc)  

  # Function to be evalulated for each pair: build distance histogram
  function potential(x,y,i,j,d2,u,mass) 
    d = sqrt(d2)
    u = u - 9.8*mass[i]*mass[j]/d
    return u
  end

  # Run pairwise computation
  u = map_pairwise((x,y,i,j,d2,u) -> potential(x,y,i,j,d2,u,mass),0.0,x,box,lc)
  return u

end

#
# In this test we compute the "gravitational force", pretending that each particle
# has a different mass. In this case, the closure is used to pass the masses and
# the force vector to the function that computes the potential.
#
function test4(N=100_000)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists
  lc = LinkedLists(N)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N ]

  # masses
  mass = rand(N)

  # forces
  forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

  # Initializing linked cells with these positions
  initcells!(x,box,lc)  

  # Function to be evalulated for each pair: build distance histogram
  function calc_forces!(x,y,i,j,d2,mass,forces) 
    G = 9.8*mass[i]*mass[j]/d2
    d = sqrt(d2)
    df = (G/d) * (x - y)
    forces[i] = forces[i] - df
    forces[j] = forces[j] + df
    return forces
  end

  # Run pairwise computation
  forces = map_pairwise((x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),forces,x,box,lc)
  return forces

end

#
# In this test we compute the minimum distance between two independent sets of particles
#
function test5(;N1=1_500,N2=1_500_000)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists (largest set!)
  lc = LinkedLists(N2)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N1 ]
  y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N2 ]

  # Initializing linked cells with these positions (largest set!)
  initcells!(y,box,lc)  

  # Function that keeps the minimum distance
  f(x,y,i,j,d2,mind) = d2 < mind[3] ? mind = (i,j,d2) : mind

  # We have to define our own reduce function here
  function reduce(output_threaded::Vector{Tuple{Int,Int,Float64}})
    mind = output_threaded[1]
    for i in 2:nthreads()
      if output_threaded[i][3] < mind[3]
        mind = output_threaded[i]
      end
    end
    return (mind[1],mind[2],sqrt(mind[3]))
  end 

  # Result
  mind = ( 0, 0, +Inf )

  # Run pairwise computation
  mind = map_pairwise(f,mind,x,y,box,lc;reduce=reduce)
  return mind

end

#
# In this test we compute the minimum distance between two independent sets of particles
#
function test5(;N1=1_500,N2=1_500_000)

  # Number of particles, sides and cutoff
  sides = [250,250,250]
  cutoff = 10.
  box = Box(sides,cutoff)

  # Initialize auxiliary linked lists (largest set!)
  lc = LinkedLists(N2)

  # Particle positions
  x = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N1 ]
  y = [ box.sides .* rand(SVector{3,Float64}) for i in 1:N2 ]

  # Initializing linked cells with these positions (largest set!)
  initcells!(y,box,lc)  

  # Function that keeps the minimum distance
  f(x,y,i,j,d2,mind) = d2 < mind[3] ? mind = (i,j,d2) : mind

  # We have to define our own reduce function here
  function reduce(output_threaded::Vector{Tuple{Int,Int,Float64}})
    mind = output_threaded[1]
    for i in 2:nthreads()
      if output_threaded[i][3] < mind[3]
        mind = output_threaded[i]
      end
    end
    return (mind[1],mind[2],sqrt(mind[3]))
  end 

  # Result
  mind = ( 0, 0, +Inf )

  # Run pairwise computation
  mind = map_pairwise(f,mind,x,y,box,lc;reduce=reduce)
  return mind

end

#
# Function that uses the naive algorithm, for testing
#
function map_naive(f,output,x,box)
  cutoff2 = box.cutoff^2
  for i in 1:length(x)-1
    xᵢ = x[i]
    for j in i+1:length(x)
      xⱼ = wrapone(x[j],box.sides,xᵢ)
      d2 = sq_distance(xᵢ,xⱼ) 
      if d2 <= cutoff2
        output = f(xᵢ,xⱼ,i,j,d2,output)
      end
    end
  end
  return output
end

end # module


