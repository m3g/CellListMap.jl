"""

```
norm_sqr(v::AbstractVector{T}) where T
```

Internal function or structure - interface may change.

# Extended help

`norm_sqr` from LinearAlgebra is not documented and is slower than this for
standard arrays. 

"""
@inline function norm_sqr(v::AbstractVector{T}) where T
    n2 = zero(T)^2 # the square here is required for Units, for example
    @simd for x in v
        n2 += x^2
    end
    return n2
end

"""

```
norm(v::AbstractVector{T}) where T
```

Internal function or structure - interface may change.

# Extended help

`norm_sqr` from LinearAlgebra is not documented and is slower than this for
standard arrays. Thus we define our own `norm(x) = norm_sqr(x)`.

"""
@inline norm(v) = sqrt(norm_sqr(v))

"""

```
dot(x::AbstractVector{T1},y::AbstractVector{T2}) where {T1,T2} 
```

Internal function or structure - interface may change.

# Extended help

`LinearAlgebra.dot` is slower than this for standard arrays (likely more accurate, but
that is not relevant here).

"""
@inline function dot(x::AbstractVector{T1},y::AbstractVector{T2}) where {T1,T2} 
    length(x) == length(y) || throw(DimensionMismatch("$(length(x)) != $(length(y))"))
    d = zero(T1)*zero(T2)
    @inbounds @simd for i in eachindex(x)
        d += x[i]*y[i]
    end
    return d
end
