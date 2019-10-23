# Utilities for coding with Julia arrays

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package provides a number of methods and types to deal with the variety of
array types (sub-types of `AbstractArray`) that exist in Julia and to help
building custom array-like types.

These are useful to implement methods to process arrays in a generic way.

## Array-like objects

### Defining custom array-like objects

Julia array interface is very powerful and flexible, it is therefore tempting
to define custom array-like types, that is Julia types that behave like arrays,
without sacrificing efficiency.  The `ArrayTools` package provides simple means
to define such array-like types if the values to be accessed as if in an array
are stored in an array (of any concrete type) embedded in the object instance.

This is as simple as:

1. Make your type inherit from `LinearArray{T,N}` or `CartesianArray{T,N}`
   depending whether the index style of the embedded array is `IndexLinear()`
   or `IndexCartesian()`.

2. Extend the `Base.parent` method for your custom type so that it returns the
   embedded array.

For instance (of course replacing the ellipsis `...`):

```julia
struct CustomArray{T,N,...} <: LinearArray{T,N}
    arr::Array{T,N} # can be any array type with linear index style
    ...             # another member
    ...             # yet another member
    ...             # etc.
end

@inline Base.parent(A::CustomArray) = A.arr
```

As a result, instances of `CustomArray{T,N}` will be also seen as instances of
`AbstractArray{T,N}` and will behave as if they implement linear indexing.
Apart from the needs to extend the `Base.parent` method, the interface to
`LinearArray{T,N}` should provide any necessary methods for indexation, getting
the dimensions, the element type, *etc.* for the derived custom type.  You may
however override these definitions by more optimized or more suitable methods
specialized for your custom array-like type.

If your custom array-like type is based on an array whose index style is
`IndexCartesian()` (instead of `IndexLinear()` in the above example), just make
your custom type derived from `CartesianArray{T,N}` (instead of
`LinearArray{T,N}`).  For such array-like object, index checking requires an
efficient implementation of the `Base.axes()` method which you may have to
specialize.  The default implementation is:

```julia
@inline Base.axes(A::CartesianArray) = axes(parent(A))
```


### Array-like objects with attributes

As a working example of custom array-like objects, the `ArrayTools` package
provides `AttributeArray{T,N,K,V}` objects which store values like arrays but
also have attributes stored in a dictionary.  Here the parameters are `T` the
element type of the values in the array part, `N` the number of dimensions of
the array part, `K` the type of the attribute keys and `V` the type of the
attributes values.  To build such an object, you can do:

```julia
using ArrayTools
dims = (100, 50)
T = Float32
A = AttributeArray(zeros(T, dims), "units"=>"photons", "Δx"=>0.20, "Δy"=>0.15)
```

Then indexing `A` with integers of Cartesian indices is the same as accessing
the array of values while indexing `A` with strings is the same as accessing
the dictionary of attributes.  For example:

```julia
A["Δx"]          # yields 0.2
A["units"]       # yields "photons"
A[:,3] = 3.14    # set some values in the array part of A
sum(A)           # yeilds the sum of the values of A
A["gizmo"] = π   # set an attribute
pop!(A, "gizmo") # yields attribute value and delete it
```

An `AttributeArray` can be build by providing an array and a dictionary, an
array and key-value pairs (as in the above example), or more like other arrays:

```julia
dims = (100, 50)
T = Float32
A = AttributeArray{T}(undef, dims, "units"=>"photons", "Δx"=>0.20, "Δy"=>0.15)
```

Attribute key types are not limited to `String`, but, to avoid ambiguities,
key types must be more specialized than `Any` and must not inherit
from `Integer` nor `CartesianIndex`.

Similar types are provided by
[MetaArrays](https://github.com/haberdashPI/MetaArrays.jl),
[MetadataArrays](https://github.com/piever/MetadataArrays.jl) and
[ImageMetadata](https://github.com/JuliaImages/ImageMetadata.jl).


## General tools

### Array indexing

The `eachindex` method is very useful when writing loops over array elements so
as to be agnostic to which specfic indexing rule is the most suitable.  Some
algorithms are however more efficient or easier to write if all involved arrays
are indexed by a single 1-based index.  In that case, `using ArrayTools` provides:

```julia
fastarray(A)
```

which checks whether array `A` is suitable for fast indexing (by a single
integer starting at 1); if it does, `A` is returned to the caller; otherwise,
the contents of `A` is converted to a suitable array type implementing fast
indexing and is returned to the caller.

To just check whether array `A` is suitable for fast indexing, call:

```julia
isfastarray(A) -> bool
```

Multiple arguments can be checked at the same time: `isfastarray(A,B,...)` is
the same as `isfastarray(A) && isfastarray(B) && isflatarray(...)`.


### Array storage

When calling (with `ccall`) a compiled function coded in another language (C,
FORTRAN, etc.), you have to make sure that array arguments have the same
storage as assumed by these languages so that it is safe to pass the pointer of
the array to the compiled function.

Typically, you want to ensure that the elements are stored in memory
contiguously and in column-major order.  This can be ckecked by calling:

```julia
isflatarray(A) -> bool
```

or, with several arguments:

```julia
isflatarray(A, B, C, ...)  -> bool
```

In order to get an array with such *flat* storage and possibly with a given
element type `T`, call:

```julia
flatarray([T = eltype(A),] A)
```

which just returns `A` if the requirements hold or converts `A` to a suitable
array form.


## FAQ

* What is the difference between `IndexStyle` (defined in base Julia) and
  `IndexingTrait` (defined in `ArrayTools`)?

  If `IndexStyle(A) === IndexLinear()`, then array `A` can be efficiently
  indexed by one integer (even if `A` is multidimensional) and column-major
  ordering is used to access the elements of `A`.  The only (known) other
  possibility is `IndexStyle(A) === IndexCartesian()`.

  If `IndexingTrait(A) === FastIndexing()`, then `IndexStyle(A) ===
  IndexLinear()` also holds (see above) **and** array `A` has standard 1-based
  indices.

* What is the difference between `Base.has_offset_axes` (provided by Julia) and
  `has_standard_indexing` (provided by  `ArrayTools`)?

  For the caller, `has_standard_indexing(args...)` yields the opposite result
  as `Base.has_offset_axes(args...)`.  Furthermore, `has_standard_indexing` is
  a bit faster.

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/ArrayTools.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/ArrayTools.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/ArrayTools.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/ArrayTools.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/ArrayTools-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/ArrayTools.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/ArrayTools.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/ArrayTools.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/ArrayTools.jl?branch=master
