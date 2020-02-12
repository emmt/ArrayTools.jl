#
# traits.jl --
#
# Trait-like methods and types.
#

"""
```julia
StorageType(A)
```

yields the type of storage of the elements of argument `A`.  If `A` is a *flat*
array, that is an array with contiguous elements in column-major order and
first element at index 1, the singleton `FlatStorage()` is returned; otherwise,
the singleton `AnyStorage()` is returned.

This method can be extended for custom array types to quickly return the
correct answer.

See also [`isflatarray`](@ref), [`flatarray`](@ref).

"""
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
```julia
isflatarray(A) -> boolean
```

yields whether array `A` can be indexed as a *flat* array, that is an array
with contiguous elements in column-major order and first element at index 1.
This also means that `A` has 1-based indices along all its dimensions.

Several arguments can be checked in a single call:

```julia
isflatarray(A, B, C, ...)
```

is the same as:

```julia
isflatarray(A) && isflatarray(B) && isflatarray(C) && ...
```

See also [`StorageType`](@ref), [`flatarray`](@ref), [`isfastarray`](@ref),
[`has_standard_indexing`](@ref).

"""
isflatarray() = false
isflatarray(args...) = allof(isflatarray, args...)
isflatarray(arg) = _isflatarray(StorageType(arg))
_isflatarray(::FlatStorage) = true
_isflatarray(::StorageType) = false
#
# `isflatarray(args...)` could have been:
#
#     isflatarray(args...) = all(isflatarray, args)
#
# but using `all(isflatarray, A, x, y)` for `A`, `x` and `y` flat arrays of
# sizes (3,4,5,6), (5,6) and (3,4) takes 9.0ns with Julia 1.0 (29.1ns with
# Julia 0.6) while using `allof` takes 0.02ns (i.e. is eliminated by the
# compiler).
#

"""
```julia
flatarray([T=eltype(A),] A)
```

lazily yields a *flat* array based on `A`, that is an array with contiguous
elements in column-major order and first element at index 1.  Optional argument
`T` is to specify the element type of the result.  Argument `A` is returned if
it is already a flat array with the requested element type; otherwise,
[`convert`](@ref) is called to produce the result (an `Array{T}` in that case).

See also [`isflatarray`](@ref), [`fastarray`](@ref), [`convert`](@ref).
"""
flatarray(A::Array) = A
flatarray(::Type{T}, A::Array{T,N}) where {T,N} = A
flatarray(A::AbstractArray{T,N}) where {T,N} =
    _flatarray(StorageType(A), T, A)
flatarray(::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    _flatarray(StorageType(A), T, A)

_flatarray(::FlatStorage, ::Type{T}, A::AbstractArray{T,N}) where {T,N} = A
_flatarray(::StorageType, ::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    convert(Array{T,N}, A)

"""

```julia
IndexingTrait(A)
```

yields one of the singletons `FastIndexing()` or `AnyIndexing()` to indicate
whether or not array `A` has standard 1-based indices and can be efficiently
indexed by one integer (even if `A` is multidimensional) and column-major
ordering is used to access the elements of `A`.

This method can be extended for custom array types to quickly return the
correct answer.

See also [`isfastarray`](@ref), [`fastarray`](@ref).

"""
abstract type IndexingTrait end

struct FastIndexing <: IndexingTrait end
struct  AnyIndexing <: IndexingTrait end

@doc @doc(IndexingTrait) FastIndexing
@doc @doc(IndexingTrait) AnyIndexing

IndexingTrait(::Array) = FastIndexing()
IndexingTrait(::Any) = AnyIndexing()
IndexingTrait(A::AbstractArray) =
    (IndexStyle(A) === IndexLinear() && has_standard_indexing(A) ?
     FastIndexing() : AnyIndexing())

"""

```julia
isfastarray(A)
```

yields whether array `A` has standard 1-based indices and is efficiently
indexed by linear indices.

Several arguments can be checked in a single call:

```julia
isfastarray(A, B, C, ...)
```

is the same as:

```julia
isfastarray(A) && isfastarray(B) && isfastarray(C) && ...
```

See also [`IndexingType`](@ref), [`fastarray`](@ref), [`isflatarray`](@ref).

"""
isfastarray() = false
isfastarray(args...) = allof(isfastarray, args...)
isfastarray(arg) = _isfastarray(IndexingTrait(arg))
_isfastarray(::FastIndexing) = true
_isfastarray(::AnyIndexing) = false

"""

```julia
fastarray([T=eltype(A),] A)
```

lazily yields a *fast array* equivalent to `A` with element type `T`.  A *fast
array* has standard 1-based indices and is efficiently indexed by linear
indices.  If `A` is already a *fast array* with element type `T`, `A` is
returned; otherwise, `A` is converted into an `Array` which is returned.

See also [`isfastarray`](@ref), [`IndexingTrait`](@ref), [`flatarray`](@ref).

"""
fastarray(A::AbstractArray{T,N}) where {T,N} = fastarray(T, A)
fastarray(::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    _fastarray(IndexingTrait(A), T, A)
_fastarray(::FastIndexing, ::Type{T}, A::AbstractArray{T,N}) where {T,N} = A
_fastarray(::IndexingTrait, ::Type{T}, A::AbstractArray{<:Any,N}) where {T,N} =
    convert(Array{T,N}, A)
