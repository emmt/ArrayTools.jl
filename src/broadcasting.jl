#
# broadcasting.jl --
#
# Broadcasting of arrays with optional element type conversion.
#

"""

```julia
bcastcopy(x, [T=eltype(x),] dims...)
```

yields a new array of element type `T` and dimensions `dims` filled with `x`
according to type conversion and broadcasting rules (see [`broadcast`](@ref)).
Compared to [`bcastlazy`](@ref), it is guaranteed that the returned array does
not share its contents with `x`.

Argument `x` can be a scalar value or an array.

See also [`bcastlazy`](@ref), [`bcastdims`](@ref), [`reshape`](@ref).

"""
function bcastcopy(x, ::Type{T}, dims::Tuple{Vararg{Int}}) where {T}
    A = Array{T}(undef, dims)
    A .= x # This expression will clash if dimensions are not compatible.
    return A
end

# Convert dimensions.
bcastcopy(x, ::Type{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastcopy(x, T, map(Int, dims))
bcastcopy(x, ::Type{T}, dims::Integer...) where {T} =
    bcastcopy(x, T, dims)

# Guess element type.
bcastcopy(x, dims::Tuple{Vararg{Integer}}) =
    bcastcopy(x, eltype(x), dims)
bcastcopy(x, dims::Integer...) =
    bcastcopy(x, dims)

"""

```julia
bcastlazy(x, [T=eltype(x),] dims...)
```

yields a *flat* array of type `T` and dimensions `dims` filled with `x`
according to type conversion and broadcasting rules (see [`broadcast`](@ref)).
Compared to [`bcastcopy`](@ref), making a copy of `x` is avoided if it is
already an array with the correct type of elements and dimensions or if it can
be reshaped (see [`reshape`](@ref)) to the correct type and dimensions.  This
means that the result may share the same contents as `x`.  Argument `x` can be
a scalar or an array with 1-based indices.  The result has 1-based indices and
contiguous elements which is suitable for fast linear indexing.

See also [`bcastcopy`](@ref), [`bcastdims`](@ref), [`reshape`](@ref).

"""
function bcastlazy(A::DenseArray{T},
                   ::Type{T},
                   dims::Tuple{Vararg{Int}}) where {T}
    has_standard_indexing(A) || throw_non_standard_indexing()
    Adims = size(A)
    Adims == dims && return A
    bcastdims(Adims, dims) == dims ||
        throw(DimensionMismatch("array has incompatible dimensions"))
    return (length(A) == prod(dims) ? reshape(A, dims) : bcastcopy(A, T, dims))
end

# Call `bcastcopy` if argument `x` can certainly not be returned.
bcastlazy(x, ::Type{T}, dims::Tuple{Vararg{Int}}) where {T} =
    bcastcopy(x, T, dims)

# Convert dimensions.
bcastlazy(x, ::Type{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastlazy(x, T, map(Int, dims))
bcastlazy(x, ::Type{T}, dims::Integer...) where {T} =
    bcastlazy(x, T, dims)

# Guess element type.
bcastlazy(x, dims::Tuple{Vararg{Integer}}) =
    bcastlazy(x, eltype(x), dims)
bcastlazy(x, dims::Integer...) =
    bcastlazy(x, dims)

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
