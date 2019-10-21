# Utilities for coding with Julia arrays

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package provides a number of methods and types to deal with the variety of
array types (sub-types of `AbstractArray`) that exist in Julia.

These are useful to implement methods to process arrays in a generic way.

## Provided tools

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
