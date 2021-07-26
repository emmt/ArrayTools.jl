# Array storage

When calling (with `ccall`) a compiled function coded in another language (C,
FORTRAN, etc.), you have to make sure that array arguments have the same
storage as assumed by these languages so that it is safe to pass the pointer of
the array to the compiled function.

To check whether the elements of an array `A` are stored in memory
contiguously and in column-major order, call:

```julia
is_flat_array(A) -> bool
```

which yield a bollean result.  Several arguments can be checked in a single call:

```julia
is_flat_array(A, B, C, ...)
```

is the same as:

```julia
is_flat_array(A) && is_fast_array(B) && is_fast_array(C) && ...
```

In order to get an array with such *flat* storage and possibly with a given
element type `T`, call the [`to_flat_array`](@ref) method:

```julia
to_flat_array([T = eltype(A),] A)
```

which just returns `A` if the requirements hold for it, or which converts `A`
to a suitable array form.
