@testmodule AllocTest begin
    # This module defines the Allocs struct and the comparison operators
    # to conditionally compare the number of allocations based on the
    # BUILD_IS_PRODUCTION_BUILD environment variable.
    export Allocs
    @kwdef struct Allocs
        prodbuild::Bool = haskey(ENV, "BUILD_IS_PRODUCTION_BUILD") && ENV["BUILD_IS_PRODUCTION_BUILD"] == "true"
        allocs::Int
    end
    Allocs(allocs::Int) = Allocs(; allocs)
    import Base: ==, >, <
    ==(a::Int, b::Allocs) = b.prodbuild ? a == b.allocs : true
    <(a::Int, b::Allocs) = b.prodbuild ? a < b.allocs : true
    ==(a::Allocs, b::Int) = a.prodbuild ? a.allocs == b : true
    <(a::Allocs, b::Int) = a.prodbuild ? a.allocs < b : true
end
