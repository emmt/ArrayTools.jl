#
# traits.jl --
#
# Trait-like methods and types.
#

"""

    StorageType(A)

yields the type of storage of the elements of argument `A`.  If `A` is a *flat*
array, that is an array with contiguous elements in column-major order and
first element at index 1, the singleton `FlatStorage()` is returned; otherwise,
the singleton `AnyStorage()` is returned.

This method can be extended for custom array types to quickly return the
correct answer.

See also [`is_flat_array`](@ref), [`to_flat_array`](@ref).

""" StorageType

abstract type StorageType end

struct FlatStorage <: StorageType end
struct  AnyStorage <: StorageType end

@doc @doc(StorageType) FlatStorage
@doc @doc(StorageType)  AnyStorage

StorageType(::Array) = FlatStorage()
StorageType(::Any) = AnyStorage()
StorageType() = AnyStorage()

StorageType(A::StridedVector) =
    (first(axes(A,1)) == 1 && stride(A,1) == 1) ? FlatStorage() : AnyStorage()
StorageType(A::StridedMatrix) = begin
    inds, dims, stds = axes(A), size(A), strides(A)
    (first(inds[1]) == 1 && stds[1] == 1 &&
     first(inds[2]) == 1 && stds[2] == dims[1]) ?
     FlatStorage() : AnyStorage()
end
StorageType(A::StridedArray{T,3}) where {T} = begin
    inds, dims, stds = axes(A), size(A), strides(A)
    (first(inds[1]) == 1 && stds[1] == 1 &&
     first(inds[2]) == 1 && stds[2] == dims[1] &&
     first(inds[3]) == 1 && stds[3] == dims[1]*dims[2]) ?
     FlatStorage() : AnyStorage()
end
StorageType(A::StridedArray{T,4}) where {T} = begin
    inds, dims, stds = axes(A), size(A), strides(A)
    (first(inds[1]) == 1 && stds[1] == 1 &&
     first(inds[2]) == 1 && stds[2] == dims[1] &&
     first(inds[3]) == 1 && stds[3] == dims[1]*dims[2] &&
     first(inds[4]) == 1 && stds[4] == dims[1]*dims[2]*dims[3]) ?
     FlatStorage() : AnyStorage()
end
StorageType(A::StridedArray{T,N}) where {T,N} = begin
    inds, dims, stds = axes(A), size(A), strides(A)
    n = 1
    @inbounds for d in 1:N
        if first(inds[d]) != 1 || stds[d] != n
            return AnyStorage()
        end
        n *= dims[d]
    end
    return FlatStorage()
end

"""

    is_flat_array(A) -> boolean

yields whether array `A` can be indexed as a *flat* array, that is an array
with contiguous elements in column-major order and first element at index 1.
This also means that `A` has 1-based indices along all its dimensions.

Several arguments can be checked in a single call:

    is_flat_array(A, B, C, ...)

is the same as:

```julia
is_flat_array(A) && is_flat_array(B) && is_flat_array(C) && ...
```

See also [`StorageType`](@ref), [`to_flat_array`](@ref),
[`is_fast_array`](@ref), [`has_standard_indexing`](@ref).

"""
is_flat_array() = false
is_flat_array(args...) = allof(is_flat_array, args...)
is_flat_array(arg) = _is_flat_array(StorageType(arg))
_is_flat_array(::FlatStorage) = true
_is_flat_array(::StorageType) = false
#
# `is_flat_array(args...)` could have been:
#
#     is_flat_array(args...) = all(is_flat_array, args)
#
# but using `all(is_flat_array, A, x, y)` for `A`, `x` and `y` flat arrays of
# sizes (3,4,5,6), (5,6) and (3,4) takes 9.0ns with Julia 1.0 (29.1ns with
# Julia 0.6) while using `allof` takes 0.02ns (i.e. is eliminated by the
# compiler).
#

"""

    to_flat_array([T=eltype(A),] A)

lazily yields a *flat* array based on `A`, that is an array with contiguous
elements in column-major order and first element at index 1.  Optional argument
`T` is to specify the element type of the result.  Argument `A` is returned if
it is already a flat array with the requested element type; otherwise,
`convert` method is called to produce the result (an `Array{T}` in that case).

See also [`is_flat_array`](@ref), [`to_fast_array`](@ref).

"""
to_flat_array(A::Array) = A
to_flat_array(::Type{T}, A::Array{T,N}) where {T,N} = A
to_flat_array(A::AbstractArray{T,N}) where {T,N} =
    _to_flat_array(StorageType(A), T, A)
to_flat_array(::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    _to_flat_array(StorageType(A), T, A)

_to_flat_array(::FlatStorage, ::Type{T}, A::AbstractArray{T,N}) where {T,N} = A
_to_flat_array(::StorageType, ::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    convert(Array{T,N}, A)

"""

    IndexingType(A)

yields one of the singletons `FastIndexing()` or `AnyIndexing()` to indicate
whether or not array `A` has standard 1-based indices and can be efficiently
indexed by one integer (even if `A` is multidimensional) and column-major
ordering is used to access the elements of `A`.

This method can be extended for custom array types to quickly return the
correct answer.

See also [`is_fast_array`](@ref), [`to_fast_array`](@ref).

""" IndexingType

abstract type IndexingType end

struct FastIndexing <: IndexingType end
struct  AnyIndexing <: IndexingType end

@doc @doc(IndexingType) FastIndexing
@doc @doc(IndexingType) AnyIndexing

IndexingType(::Array) = FastIndexing()
IndexingType(::Any) = AnyIndexing()
IndexingType(A::AbstractArray) =
    (IndexStyle(A) === IndexLinear() && has_standard_indexing(A) ?
     FastIndexing() : AnyIndexing())

"""

    is_fast_array(A)

yields whether array `A` has standard 1-based indices and is efficiently
indexed by linear indices.

Several arguments can be checked in a single call:

    is_fast_array(A, B, C, ...)

is the same as:

    is_fast_array(A) && is_fast_array(B) && is_fast_array(C) && ...

See also [`IndexingType`](@ref), [`to_fast_array`](@ref),
[`is_flat_array`](@ref).

"""
is_fast_array() = false
is_fast_array(args...) = allof(is_fast_array, args...)
is_fast_array(arg) = _is_fast_array(IndexingType(arg))
_is_fast_array(::FastIndexing) = true
_is_fast_array(::AnyIndexing) = false

"""

    to_fast_array([T=eltype(A),] A)

lazily yields a *fast array* equivalent to `A` with element type `T`.  A *fast
array* has standard 1-based indices and is efficiently indexed by linear
indices.  If `A` is already a *fast array* with element type `T`, `A` is
returned; otherwise, `A` is converted into an `Array` which is returned.

See also [`is_fast_array`](@ref), [`IndexingType`](@ref),
[`to_flat_array`](@ref).

"""
to_fast_array(A::AbstractArray{T,N}) where {T,N} = to_fast_array(T, A)
to_fast_array(::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    _to_fast_array(IndexingType(A), T, A)
_to_fast_array(::FastIndexing, ::Type{T}, A::AbstractArray{T,N}) where {T,N} = A
_to_fast_array(::IndexingType, ::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    convert(Array{T,N}, A)
