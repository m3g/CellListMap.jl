#
# Some internal function to help the `show` functions to be prettier
#

# Adapt round to adhere with units
_uround(x) = round(x / oneunit(x); digits = 2) * oneunit(x)

# Get identiation from parent struct in printing
# see: https://discourse.julialang.org/t/show-nested-struct/85298/3?u=lmiq
_print(io, args...) = print(io, ' '^get(io, :indent, 0), args...)
_println(io, args...) = println(io, ' '^get(io, :indent, 0), args...)
