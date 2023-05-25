# Broadcasting of arrays

`ArrayTools` provides the [`bcastlazy`](@ref) and [`bcastcopy`](@ref) methods
for broadcasting arrays (that is expand an array to change the lengths of unit
dimensions or to add trailing dimensions) with optional element type
conversion.

For instance:

```julia
bcastlazy(A, [T=eltype(A),] dims...)
```

yields a *flat* array (see [`is_flat_array`](@ref)) of type `T` and dimensions
`dims` whose values are given by `A` and according to type conversion and
broadcasting rules (see `broadcast` standard method).

Making a copy of `A` is avoided by [`bcastlazy`](@ref) (hence its name), and
`A` or a view on `A` is returned if `A` is already an array with the correct
element type and dimensions or if it can be reshaped (by the `reshape` method)
to match the constraints. This means that the result may share the same
contents as `A`. If an array that does not share its contents with `A` is
needed, then call:

```julia
bcastcopy(A, [T=eltype(A),] dims...)
```

instead of  [`bcastlazy`](@ref).

The result returned by [`bcastlazy`](@ref) and [`bcastcopy`](@ref) has 1-based
indices and contiguous elements which is suitable for fast linear indexing.

The [`bcastsize`](@ref) method may be useful to determine dimensions when
applying broadcasting rules:

```julia
bcastsize(size(A), size(B), ...) -> siz
bcastsize(a, b) -> c
```

The first one yields the size `siz` of the array that would result from
applying broadcasting rules to arguments `A`, `B`, etc. The result is a tuple
of integers (of type `Int`). You may call [`check_size`](@ref) to ensure that
the result is a valid list on nonnegative dimensions.

The second applies broadcasting rules for a single dimension, throwing an
exception if dimensions `a` and `b` are not compatible according to these
rules.
