# Utilities for coding with Julia arrays

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This [Julia][julia-url] package provides a number of methods and types to deal
with the variety of array types (sub-types of `AbstractArray`) that exist in
Julia and to help building custom array-like types without sacrificing
performances.

These are useful to implement methods to process arrays in a generic way.

## Rubber indices

The constant `…` (type `\ldots` and hit the `tab` key) can be used in array
indexation to left or right justify the other indices.  For instance, assuming
`A` is a `3×4×5×6` array, then all the following equalities hold:

```julia
A[…]           == A[:,:,:,:]
A[…,3]         == A[:,:,:,3]
A[2,…]         == A[2,:,:,:]
A[…,2:4,5]     == A[:,:,2:4,5]
A[:,2:3,…]     == A[:,2:3,:,:]
A[2:3,…,1,2:4] == A[2:3,:,1,2:4]
```

As you can see, the advantage of the *rubber index* `…` is that it
automatically expands as the number of colons needed to have the correct number
of indices.  The expressions are also more readable.  The idea comes from the
[`Yorick`](http://yorick.github.com/) language by Dave Munro.

The rubber index may also be used for setting values.  For instance:

```julia
A[…] .= 1        # to fill A with ones
A[…,3] = A[…,2]  # to copy A[:,:,:,2] in A[:,:,:,3]
A[…,3] .= A[…,2] # idem but faster
A[2,…] = A[3,…]  # to copy A[3,:,:,:] in A[2,:,:,:]
A[…,2:4,5] .= 7  # to set all elements in A[:,:,2:4,5] to 7
```

Leading/trailing indices may be specified as Cartesian indices (of type
`CartesianIndex`).

Technically, the constant `…` is defined as `RubberIndex()` where `RubberIndex`
is the singleron type that represents any number of indices.

Call `rubberindex(n)` if you need a *rubber index* of length `n`, that is a
`n`-tuple of colons. When `n` is known at compile time, it is faster to call
`rubberindex(Val(n))`.


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

2. Extend the `Base.parent(A)` method for your custom type so that it returns
   the embedded array of an instance `A`.

For instance (of course replacing the ellipsis `...`):

```julia
using ArrayTools.PseudoArrays
struct CustomArray{T,N,...} <: LinearArray{T,N}
    arr::Array{T,N} # can be any array type with linear index style
    ...             # another member
    ...             # yet another member
    ...             # etc.
end

@inline Base.parent(A::CustomArray) = A.arr
```

As a result, instances of your `CustomArray{T,N}` will be also seen as
instances of `AbstractArray{T,N}` and will behave as if they implement linear
indexing.  Apart from the needs to extend the `Base.parent` method, the
interface to `LinearArray{T,N}` should provide any necessary methods for
indexation, getting the dimensions, the element type, *etc.* for the derived
custom type.  You may however override these definitions by more optimized or
more suitable methods specialized for your custom array-like type.

If your custom array-like type is based on an array whose index style is
`IndexCartesian()` (instead of `IndexLinear()` in the above example), just make
your custom type derived from `CartesianArray{T,N}` (instead of
`LinearArray{T,N}`).  For such array-like object, index checking requires an
efficient implementation of the `Base.axes()` method which you may have to
specialize.  The default implementation is:

```julia
@inline Base.axes(A::CartesianArray) = axes(parent(A))
```


### Array-like objects with properties

As a working example of custom array-like objects, the `ArrayTools` package
provides `AnnotatedArray{T,N,P}` objects which store values like arrays but
also have properties stored in a dictionary or a named tuple (of type `P`).
Here the parameters are the element type `T` of the values in the array part,
the number `N` of dimensions of the array part and the type `P` of the object
storing the properties.

Building annotated arrays is easy:

```julia
using ArrayTools.AnnotatedArrays
dims = (100, 50)
T = Float32
A = AnnotatedArray(zeros(T, dims), units = "photons", Δx = 0.20, Δy = 0.15)
B = AnnotatedArray{T}(undef, dims, units = "µm", Δx = 0.10, Δy = 0.20)
```

Here the initial properties of `A` and `B` are specified by the keywords in the
call to the constructor; their properties will have symbolic names with any
kind of value.  The array contents of `A` is an array of zeros, while the array
contents of `B` is created by the constructor with undefined values.  Indexing
`A` or `B` with integers of Cartesian indices is the same as accessing the
values of their array contents while indexing `A` or `B` by symbols is the same
as accessing their properties.  For example:

```julia
A.Δx             # yields 0.2
A[:Δx]           # idem
A.units          # yields "photons"
A[:units]        # idem
A[:,3] .= 3.14   # set some values in the array contents of A
sum(A)           # yields the sum of the values of A
A[:gizmo] = π    # set a property
A.gizmo = π      # idem
pop!(A, :gizmo)  # yields property value and delete it
```

Creating an annotated array is summarized by:

```julia
using ArrayTools.AnnotatedArrays
A = AnnotatedArray(arr, prop)
B = AnnotatedArray{T}(init, dims, prop)
```

where `arr` is an existing array or an expression whose result is an array,
`prop` specifies the initial properties (more on this below), `T` is the type
of array element, `init` is usually `undef` and `dims` is a tuple of array
dimensions.  If `arr` is an existing array, the object `A` created above will
reference this array and hence share its contents with the caller (call
`copy(arr)` to avoid that).  The same applies if the initial properties are
specified by a dictionary.

The properties `prop` can be specified by keywords, by key-value pairs, as a
dictionary or as a named tuple.  To avoid ambiguities, these different styles
cannot be mixed.  Below are a few examples:

```julia
using ArrayTools.AnnotatedArrays
arr = zeros(3,4,5)
A = AnnotatedArray(arr,  units  =  "µm",  Δx  =  0.1,  Δy  =  0.2)
B = AnnotatedArray(arr, :units  => "µm", :Δx  => 0.1, :Δy  => 0.2)
C = AnnotatedArray(arr, "units" => "µm", "Δx" => 0.1, "Δy" => 0.2)
D = AnnotatedArray(arr, (units  =  "µm",  Δx  =  0.1,  Δy  =  0.2))
```

The two first examples (`A` and `B`) both yield an annotated array whose
properties have symbolic keys and can have any type of value.  The third
example (`C`) yields an annotated array whose properties have string keys and
can have any type of value.  The properties of `A`, `B` and `C` are *dynamic*:
they can be modified, deleted and new properties can be inserted.  The fourth
example (`D`) yields an annotated array whose properties are stored by a *named
tuple*, they are *immutable* and have symbolic keys.

Accessing a property is possible via the syntax `obj[key]` or, for symbolic and
textual keys, via the syntax `obj.key`.  Accessing *immutable* properties is
the fastest while accessing textual properties as `obj.key` is the slowest
(because it involves converting a symbol into a string).

When initially specified by keywords or as key-value pairs, the properties are
stored in a dictionary whose key type is specialized if possible (for
efficiency) but with value type `Any` (for flexibility).  If one wants specific
properties key and value types, it is always possible to explicitly specify a
dictionary in the call to `AnnotatedArray`.  For instance:

```julia
E = AnnotatedArray(arr, Dict{Symbol,Int32}(:a => 1, :b => 2))
```

yields an annotated array whose properties have symbolic keys and integer
values of type `Int32`.

Property key types are not limited to `Symbol` or `String`, but, to avoid
ambiguities, key types must be more specialized than `Any` and must not inherit
from types like `Integer` or `CartesianIndex` which are reserved for indexing
the array contents of annotated arrays.

If the dictionary is unspecified, the properties are stored in a, initially
empty, dictionary with symbolic keys and value of any type, *i.e.*
`Dict{Symbol,Any}()`.

Iterating on an annotated array is iterating on its array values.  To iterate
on its properties, call the `properties` method which returns the object
storing the properties:

```julia
dims = (100, 50)
T = Float32
N = length(dims)
A = AnnotatedArray(zeros(T, dims), units = "µm", Δx = 0.2, Δy = 0.1)
for (k,v) in properties(A)
    println(k, " => ", v)
end
```

Similar types are provided by
[MetaArrays](https://github.com/haberdashPI/MetaArrays.jl),
[MetadataArrays](https://github.com/piever/MetadataArrays.jl) and
[ImageMetadata](https://github.com/JuliaImages/ImageMetadata.jl).


## General tools

### Array indexing

The `safeindices` method takes any number of array arguments and yields an
efficient iterator for visiting all indices each index of the arguments.  Its
behavior is similar to that of `eachindex` method except that `safeindices`
throws a `DimensionMismatch` exception if the arrays have different axes.  As a
consequence, it is always safe to specify `@inbounds` for a loop like:

```julia
@inbounds for i in safeindices(A, B, C, D)
   A[i] = B[i]*C[i] + D[i]
end
```

The `eachindex` and `safeindices` methods are very useful when writing loops
over array elements so as to be agnostic to which specfic indexing rule is the
most suitable.  Some algorithms are however more efficient or easier to write
if all involved arrays are indexed by a single 1-based index.  In that case,
`using ArrayTools` provides:

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


## Installation

`ArrayTools` is not yet an [official Julia package][julia-pkgs-url] so you have
to clone its repository.  In Julia, hit the `]` key to switch to the package
manager REPL (you should get a `... pkg>` prompt) and type:

```julia
... pkg> add https://github.com/emmt/ArrayTools.jl
```


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

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/
