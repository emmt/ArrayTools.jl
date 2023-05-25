# Array Indexing

The [`all_indices`](@ref) method takes any number of array arguments and yields
an efficient iterator for visiting all indices each index of the arguments. Its
behavior is similar to that of `eachindex` method except that `all_indices`
throws a `DimensionMismatch` exception if linearly indexed arrays have
different axes while `eachindex` just checks that such arrays have the same
number of elements. It is always safe to specify `@inbounds` (and `@simd` or
`@turbo` from the `LoopVectorization` package) for a loop like:

```julia
for i in all_indices(A, B, C, D)
   A[i] = B[i]*C[i] + D[i]
end
```

An alternative is to call the [`@assert_same_axes`](@ref) macro which throws a
`DimensionMismatch` exception if the provided arguments are not arrays with the
same axes. For example:

```julia
@assert_same_axes A B C D
@inbounds for i in eachindex(A, B, C, D)
   A[i] = B[i]*C[i] + D[i]
end
```

where the macro call amounts to:

```julia
axes(A) == axes(B) == axes(C) == axes(D) ? nothing : throw(DimensionMismatch("..."))
```

The `eachindex` and [`all_indices`](@ref) methods are very useful when writing
loops over array elements so as to be agnostic to which specific indexing rule
is the most suitable. Some algorithms are however more efficient or easier to
write if all involved arrays are indexed by a single 1-based index.
`ArrayTools` provides [`is_fast_array`](@ref) to check whether arrays are
suitable for fast indexing:

```julia
is_fast_array(A) -> bool
```

to check array `A` or:

```julia
is_fast_array(A, B, ...) -> bool
```

to check several arrays `A`, `B`, ... at the same time. `ArrayTools` also
provides [`to_fast_array`](@ref) to convert an array to fast linear indexing if
this is needed:

```julia
to_fast_array(A)
```

checks whether array `A` is suitable for fast indexing (by a single integer
starting at 1); if it does, `A` is returned to the caller; otherwise, the
contents of `A` is converted to a suitable array type implementing fast
indexing and is returned to the caller.

The method [`has_standard_indexing`](@ref) can be called to check whether one
or several arrays have standard 1-based indexing though they not necessarily
implement fast linear indexing.
