module ArrayTools

export
    …,
    RubberIndex,
    allof,
    anyof,
    bcastcopy,
    bcastdim,
    bcastdims,
    bcastlazy,
    cartesian_indices,
    checkdimensions,
    colons,
    Dimensions,
    dimensions,
    has_standard_indexing,
    noneof,
    promote_eltype,
    reversemap,
    safe_indices,
    # storage trait
    StorageType,
    AnyStorage,
    FlatStorage,
    flatarray,
    isflatarray,
    # Fast arrays and indexing trait:
    IndexingTrait,
    FastIndexing,
    AnyIndexing,
    fastarray,
    isfastarray

using Base: OneTo, axes1, @_inline_meta
import Base: dotview, getindex, setindex!

@deprecate rubberindex colons
@deprecate indices cartesian_indices
@deprecate cartesianindices cartesian_indices
@deprecate safeindices safe_indices
@deprecate flatmatrix flatarray
@deprecate flatvector flatarray
@deprecate fastmatrix fastarray
@deprecate fastvector fastarray

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

"""

```julia
promote_eltype(A, B, ...) -> T
```

yields an element type `T` resulting from calling `promote_type` onto the
element types of all the arguments which must be arrays or array types.

"""
promote_eltype(args::Union{AbstractArray,Type{<:AbstractArray}}...) =
    promote_type(map(eltype, args)...)

"""

```julia
dimensions(arg) -> dims
```

yields a list of dimensions (as a tuple of `Int`) out of argument `arg`.  The
union `Dimensions` matches the types of argument `arg` acceptable for
`dimensions(arg)`: scalar integer and tuple of integers.

As a special case, for an array `A`:

```julia
dimensions(A) -> size(A)
```

yields the list of dimensions of `A`, that is `size(A)`, after having checked
that `A` has standard 1-based indices.

See also [`has_standard_indexing`](@ref).

"""
dimensions(dims::Tuple{}) = dims
dimensions(dim::Integer) = (Int(dim),)
dimensions(dims::Tuple{Vararg{Integer}}) = map(Int, dims)
dimensions(dims::Integer...) = map(Int, dims)
dimensions(dim::Int) = (dim,)
dimensions(dims::Tuple{Vararg{Int}}) = dims
dimensions(dims::Int...) = dims

function dimensions(A::AbstractArray)
    has_standard_indexing(A) || throw_non_standard_indexing()
    return size(A)
end

const Dimensions = Union{Integer,Tuple{Vararg{Integer}}}

@doc @doc(dimensions) Dimensions

@noinline throw_non_standard_indexing() =
    error("array have non-standard indices")

"""

```julia
checkdimensions(dims)
```

checks whether `dims` is a list of valid dimensions.  An error is thrown if not
all dimensions are nonnegative.

See also [`dimensions`](@ref).

"""
checkdimensions(::Tuple{}) = true
checkdimensions(dims::Tuple{Vararg{Integer}}) =
    allof(d -> d ≥ 0, dims...) || error("invalid array dimension(s)")
checkdimensions(dim::Integer) =
    dim ≥ 0 || error("invalid array dimension")

"""

```julia
has_standard_indexing(A)
```

return `true` if the indices of `A` start with 1 along all axes.  Can be called
with multiple arguments:

```julia
has_standard_indexing(A, B, ...)
```

is equivalent to:

```julia
has_standard_indexing(A) && has_standard_indexing(B) && ...
```

Opposite of `Base.has_offset_axes` which is not available in version of Julia
older than 0.7.

"""
has_standard_indexing(arg) = allof(x -> first(x) == 1, axes(arg)...)
has_standard_indexing(args...) = allof(has_standard_indexing, args...)
has_standard_indexing(arg::Array) = true
has_standard_indexing(::Colon) = true
has_standard_indexing(::Number) = true

"""

The calls:

```julia
cartesian_indices(A)
cartesian_indices((n1, n2, ...))
cartesian_indices((i1:j1, i2:j2, ...))
cartesian_indices(CartesianIndex(i1, i2, ...), CartesianIndex(j1, j2, ...))
cartesian_indices(R)
```

all yield an instance of `CartesianIndices` suitable for multi-dimensional
indexing of respectively: all the indices of array `A`, a multi-dimensional
array of dimensions `(n1,n2,...)`, a multi-dimensional region whose first and
last indices are `(i1,i2,...)` and `(j1,j2,...)` or a Cartesian region defined
by `R`, an instance of `CartesianIndices`.

"""
cartesian_indices(A::AbstractArray) = cartesian_indices(axes(A))
cartesian_indices(R::CartesianIndices) = R
@inline function cartesian_indices(start::CartesianIndex{N},
                                   stop::CartesianIndex{N}) where {N}
    CartesianIndices(map((i,j) -> (i == 1 ? Base.OneTo(j) : i:j),
                         start.I, stop.I))
end
cartesian_indices(dims::Tuple{Vararg{Integer}}) =
    CartesianIndices(map(dim -> Base.OneTo(dim), dims))
cartesian_indices(rngs::NTuple{N,AbstractUnitRange{<:Integer}}) where {N} =
    CartesianIndices(rngs)

# The following, would yield an `AbstractUnitRange` if a single argument is
# provided that is an integer or a range (not a tuple).  This lead to
# ambiguities so it is has been disabled.  The rules are that (i) an array
# argument (including a range) yields the indices for this array, (ii) a tuple
# of dimensions, a tuple of unit ranges or a pair of `CartesianIndices` yields
# the indices of the Cartesian region defined by these arguments.
#
#cartesian_indices(dim::Int) = Base.OneTo(dim)
#cartesian_indices(dim::Integer) = Base.OneTo(Int(dim))
#cartesian_indices(rng::AbstractUnitRange{Int}) = rng
#cartesian_indices(rng::AbstractUnitRange{<:Integer}) = convert(UnitRange{Int}, rng)


#------------------------------------------------------------------------------

"""

```julia
safe_indices(A...)
```

yields an iterable object for visiting each index of array(s) `A` in an
efficient manner. For array types that have opted into fast linear indexing
(like `Array`), this is simply the range `1:length(A)`. For other array types,
return a specialized Cartesian range to efficiently index into the array with
indices specified for every dimension.

If more than one `AbstractArray` argument are supplied, `safe_indices` will
create an iterable object that is fast for all arguments (a `UnitRange` if all
inputs have fast linear indexing, a `CartesianIndices` otherwise).  A
`DimensionMismatch` exception is thrown if the arrays have different axes so
that it is always safe to use `@inbounds` in front of a loop like:

```julia
for i in safe_indices(A, B, C, D)
   A[i] = B[i]*C[i] + D[i]
end
```

when `A`, `B` etc. are all (abstract) arrays.

This method is similar to `eachindex` except that a `DimensionMismatch`
exception is thrown if arrays have different axes.  For linearly indexed
arrays, `eachindex` only checks that they have the same linear index range
(that is the same number of elements, not the same shape).

""" safe_indices

@inline safe_indices(A::AbstractArray) =
    eachindex(IndexStyle(A), A)
@inline safe_indices(A::AbstractArray, B::AbstractArray) =
    safe_indices(IndexStyle(A, B), A, B)
@inline safe_indices(A::AbstractArray, B::AbstractArray...) =
    safe_indices(IndexStyle(A, B...), A, B...)

@inline function safe_indices(::IndexLinear, A::AbstractArray,
                              B::AbstractArray...)
    all_match_first(axes, axes(A), B...) ||
        throw_indices_mismatch(IndexLinear(), A, B...)
    return eachindex(IndexLinear(), A)
end

@inline function safe_indices(::IndexCartesian, A::AbstractArray,
                              B::AbstractArray...)
    inds = axes(A)
    all_match_first(axes, inds, B...) ||
        throw_indices_mismatch(IndexCartesian(), A, B...)
    # The following is the same as `eachindex(IndexCartesian(),A)` which yields
    # `CartesianIndices(axes(A))`.
    return CartesianIndices(inds)
end

all_match_first(f::Function, val, arg, args...) = begin
    @_inline_meta
    (val == f(arg)) & all_match_first(f, val, args...)
end
all_match_first(f::Function, val) = true

@noinline function throw_indices_mismatch(::IndexLinear,
                                          A::AbstractArray...)
    io = IOBuffer()
    write(io, "all arguments of `safe_indices` must have the same shape, got ")
    m = length(A)
    for i in 1:m
        write(io, (i == 1 ? "(" : i == m ? " and (" : ", )"))
        dims = size(A[i])
        n = length(dims)
        for j in 1:n
            j > 1 && write(io, ",")
            print(io, dims[j])
        end
        write(io, (n == 1 ? ",)" : ")"))
    end
    throw(DimensionMismatch(String(take!(io))))
end

@noinline function throw_indices_mismatch(::IndexCartesian,
                                          A::AbstractArray...)
    io = IOBuffer()
    write(io, "all arguments of `safe_indices` must have the same axes, got ")
    m = length(A)
    for i in 1:m
        write(io, (i == 1 ? "(" : i == m ? " and (" : ", )"))
        inds = axes(A[i])
        n = length(inds)
        for j in 1:n
            j > 1 && write(io, ",")
            print(io, first(inds[j]), ":", last(inds[j]))
        end
        write(io, (n == 1 ? ",)" : ")"))
    end
    throw(DimensionMismatch(String(take!(io))))
end

#------------------------------------------------------------------------------
# BROADCASTING OF ARRAYS WITH OPTIONAL ELEMENT TYPE CONVERSION.

"""

```julia
bcastcopy([T=eltype(A),] A, dims...)
```

yields a new array of type `T` and dimensions `dims` filled with the elements
of `A` according to type conversion and broadcasting rules (see
[`broadcast`](@ref)).  Compared to [`bcastlazy`](@ref), it is guaranteed that
the returned array does not share its contents with `A`.

Arguments `A` and `T` can be exchanged, that is `bcastcopy(A,T,dims)` is the
same as `bcastcopy(T,A,dims)`.

See also [`bcastlazy`](@ref), [`bcastdims`](@ref), [`reshape`](@ref).

"""
bcastcopy(A::AbstractArray{T}, dims::Integer...) where {T} =
    bcastcopy(T, A, dims)
bcastcopy(A::AbstractArray{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastcopy(T, A, dims)
bcastcopy(x::T, dims::Integer...) where {T} =
    bcastcopy(T, x, dims)
bcastcopy(x::T, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastcopy(T, x, dims)
bcastcopy(x, ::Type{T}, dims::Integer...) where {T} =
    bcastcopy(T, x, dims)
bcastcopy(x, ::Type{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastcopy(T, x, dims)
bcastcopy(::Type{T}, x, dims::Integer...) where {T} =
    bcastcopy(T, x, dims)
bcastcopy(::Type{T}, x, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastcopy(T, x, map(Int, dims))

function bcastcopy(::Type{T}, A::AbstractArray,
                   dims::Tuple{Vararg{Int}}) where {T}
    C = Array{T}(undef, dims)
    @. C = A # This exprssion will clash if dimensions are not compatible.
    return C
end

bcastcopy(::Type{T}, x::T, dims::Tuple{}) where {T} = copy(x)
bcastcopy(::Type{T}, x, dims::Tuple{}) where {T} = convert(T, x)
bcastcopy(::Type{T}, x, dims::NTuple{N,Int}) where {N,T} =
    fill!(Array{T,N}(undef, dims), x)

"""

```julia
bcastlazy([T=eltype(A),] A, dims...)
```

yields a *flat* array of type `T` and dimensions `dims` filled with the
elements of `A` according to type conversion and broadcasting rules (see
[`broadcast`](@ref)).  Compared to [`bcastcopy`](@ref), making a copy of `A` is
avoided if it already has the correct type and dimensions or can be reshaped
(see [`reshape`](@ref)) to the correct type and dimensions.  This means that
the result may share the same contents as `A`.  Array `A` must have 1-based
indices.  The result has 1-based indices and contiguous elements which is
suitable for fast linear indexing.

Arguments `A` and `T` can be exchanged, that is `bcastlazy(A,T,dims)` is the
same as `bcastlazy(T,A,dims)`.

See also [`bcastcopy`](@ref), [`bcastdims`](@ref), [`reshape`](@ref).

"""
bcastlazy(A::AbstractArray{T}, dims::Integer...) where {T} =
    bcastlazy(T, A, dims)
bcastlazy(A::AbstractArray{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastlazy(T, A, dims)
bcastlazy(x::T, dims::Integer...) where {T} =
    bcastlazy(T, x, dims)
bcastlazy(x::T, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastlazy(T, x, dims)
bcastlazy(x, ::Type{T}, dims::Integer...) where {T} =
    bcastlazy(T, x, dims)
bcastlazy(x, ::Type{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastlazy(T, x, dims)
bcastlazy(::Type{T}, x, dims::Integer...) where {T} =
    bcastlazy(T, x, dims)
bcastlazy(::Type{T}, x, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastlazy(T, x, map(Int, dims))

function bcastlazy(::Type{T}, A::AbstractArray{T},
                   dims::Tuple{Vararg{Int}}) where {T}
    has_standard_indexing(A) || throw_non_standard_indexing()
    Adims = size(A)
    Adims == dims && return A
    bcastdims(Adims, dims) == dims ||
        throw(DimensionMismatch("array has incompatible dimensions"))
    return (isa(A, DenseArray) && length(A) == prod(dims) ? reshape(A, dims) :
            bcastcopy(T, A, dims))
end

bcastlazy(::Type{T}, A::AbstractArray, dims::Tuple{Vararg{Int}}) where {T} =
    bcastcopy(T, A, dims)

bcastlazy(::Type{T}, x::T, dims::Tuple{}) where {T} = x
bcastlazy(::Type{T}, x, dims::Tuple{}) where {T} = convert(T, x)
bcastlazy(::Type{T}, x, dims::Tuple{Vararg{Int}}) where {T} =
    bcastcopy(T, x, dims)

"""

```julia
bcastdims(size(A), size(B), ...)
```

yields the dimensions of the array that would result from applying broadcasting
rules (see [`broadcast`](@ref)) to arguments `A`, `B`, etc.  The result is a
tuple of dimensions of type `Int`.  Call [`checkdimensions`](@ref) if you want
to also make sure that the result is a list of valid dimensions.

See also [`dimensions`](@ref), [`checkdimensions`](@ref), [`bcastcopy`](@ref),
[`bcastlazy`](@ref).

"""
bcastdims(::Tuple{}) = ()
bcastdims(a::Tuple{Vararg{Integer}}) = dimensions(a)
bcastdims(a::Tuple{Vararg{Integer}}, b::Tuple{Vararg{Integer}}, args...) =
    bcastdims(bcastdims(a, b), args...)

# Use a recursion to build the dimension list from two lists, if code is
# inlined (for a few number of dimensions), it should be very fast.
bcastdims(::Tuple{}, ::Tuple{}) = ()
bcastdims(::Tuple{}, b::Tuple{Vararg{Integer}}) = bcastdims(b)
bcastdims(a::Tuple{Vararg{Integer}}, ::Tuple{}) = bcastdims(a)
bcastdims(a::Tuple{Vararg{Integer}}, b::Tuple{Vararg{Integer}}) =
    (bcastdim(a[1], b[1]), bcastdims(Base.tail(a), Base.tail(b))...)

# Apply broadcasting rules for a single dimension (same as
# Base.Broadcasting._bcs1 but takes care of converting to `Int`).
bcastdim(a::Integer, b::Integer) = bcastdim(Int(a), Int(b))
bcastdim(a::Int, b::Int) =
    (a == b || b == 1 ? a : a == 1 ? b :
    throw(DimensionMismatch("arrays could not be broadcast to a common size")))

#------------------------------------------------------------------------------

"""
```julia
allof(p, args...) -> Bool
```

checks whether predicate `p` returns `true` for all arguments in `args...`,
returning `false` as soon as possible (short-circuiting).

```julia
allof(args...) -> Bool
```

checks whether all arguments `args...` are `true`, returning `false` as soon as
possible (short-circuiting).  Arguments can be booleans or arrays of booleans.
The latter are considered as `true` if all their elements are `true` and are
considered as `false` otherwise (if any of their elements are `false`).
Arguments can also be iterables to check whether all their values are `true`.
An empty iterable is considered as `true`.

This method can be much faster than `all(p, args)` or `all(args)` because its
result may be determined at compile time.  However, `missing` values are not
considered as special.

See also [`all`](@ref), [`anyof`](@ref), [`noneof`](@ref).

"""
allof(p::Function, a) = p(a)::Bool
allof(p::Function, a, b...) = p(a) && allof(p, b...)
allof(a, b...) = allof(a) && allof(b...)
allof(a::Bool) = a
function allof(a::AbstractArray{Bool})
    @inbounds for i in eachindex(a)
        a[i] || return false
    end
    return true
end
function allof(itr)
    for val in itr
        allof(val) || return false
    end
    return true
end

"""
```julia
anyof(p, args...) -> Bool
```

checks whether predicate `p` returns `true` for any argument `args...`,
returning `true` as soon as possible (short-circuiting).

```julia
anyof(args...) -> Bool
```

checks whether all arguments `args...` are `true`, returning `false` as soon as
possible (short-circuiting).  Arguments can be booleans or arrays of booleans.
The latter are considered as `true` if any of their elements are `true` and are
considered as `false` otherwise (if all their elements are `false`).  Arguments
can also be iterables to check whether any of their values are `true`.  An
empty iterable is considered as `false`.

This method can be much faster than `any(p, args)` or `any(args)` because its
result may be determined at compile time.  However, `missing` values are not
considered as special.

To check whether predicate `p` returns `false` for all argument `args...` or
whether all argument `args...` are false, repectively call:

```julia
noneof(p, args...) -> Bool
```

or

```julia
noneof(args...) -> Bool
```

which are the same as `!anyof(p, args...)` and `!anyof(args...)`.

See also [`any`](@ref), [`allof`](@ref).

"""
anyof(p::Function, a) = p(a)::Bool
anyof(p::Function, a, b...) = p(a) || anyof(p, b...)
anyof(a, b...) = anyof(a) || anyof(b...)
anyof(a::Bool) = a
function anyof(a::AbstractArray{Bool})
    @inbounds for i in eachindex(a)
        a[i] && return true
    end
    return false
end
function anyof(itr)
    for val in itr
        anyof(val) && return true
    end
    return false
end

noneof(args...) = ! anyof(args...)
@doc @doc(anyof) noneof

"""
```julia
reversemap(f, args)
```

applies the function `f` to arguments `args` in reverse order and return the
result.  For now, the arguments `args` must be in the form of a simple tuple
and the result is the tuple: `(f(args[end]),f(args[end-1]),...,f(args[1])`.

Also see: [`map`](@ref), [`ntuple`](@ref).

"""
reversemap(f::Function, args::NTuple{N,Any}) where {N} =
    ntuple(i -> f(args[(N + 1) - i]), Val(N))

include("rubberindex.jl")

include("PseudoArrays.jl")
using .PseudoArrays

include("AnnotatedArrays.jl")
using .AnnotatedArrays

end # module
