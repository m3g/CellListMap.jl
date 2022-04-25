# Help entries

These entries can be viewed at the `Julia` REPL `Julia` using 

```julia-repl
julia> ? 
help?> function_name
```

## Public interface

The types and functions provided here will have their interface unchanged among non-breaking releases.

```@autodocs
Private = false
Modules=[CellListMap]
Order = [ :type, :function ]
```

## Internals

These are internal functions, whose interface can be modified without implying a breaking release.

```@autodocs
Public = false
Modules=[CellListMap]
Order = [ :type, :function ]
```

