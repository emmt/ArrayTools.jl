#
# broadcasting.jl --
#
# Broadcasting of arrays with optional element type conversion.
#

"""

    bcastcopy(A, [T=eltype(A),] dims...)

yields a new array of element type `T` and dimensions `dims` whose values are
given by `A` according to type conversion and broadcasting rules (like for the
`broadcast` method).  Compared to [`bcastlazy`](@ref), it is guaranteed that
the returned array does not share its contents with `A`.

Argument `A` can be a scalar value or an array.

See also [`bcastlazy`](@ref), [`bcastdims`](@ref).

"""
function bcastcopy(A, ::Type{T}, dims::Tuple{Vararg{Int}}) where {T}
    B = Array{T}(undef, dims)
    B .= A # This expression will clash if dimensions are not compatible.
    return B
end

# Convert dimensions.
bcastcopy(A, ::Type{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastcopy(A, T, map(Int, dims))
bcastcopy(A, ::Type{T}, dims::Integer...) where {T} =
    bcastcopy(A, T, dims)

# Guess element type.
bcastcopy(A, dims::Tuple{Vararg{Integer}}) =
    bcastcopy(A, eltype(A), dims)
bcastcopy(A, dims::Integer...) =
    bcastcopy(A, dims)

"""

    bcastlazy(A, [T=eltype(A),] dims...)

yields a *flat* array of type `T` and dimensions `dims` whose values are given
by `A` according to type conversion and broadcasting rules (see
[`broadcast`](@ref)).  Compared to [`bcastcopy`](@ref), making a copy of `A` is
avoided if it is already an array with the correct type of elements and
dimensions or if it can be reshaped (by the `reshape` method) to the correct
type and dimensions.  This means that the result may share the same contents as
`A`.  Argument `A` can be a scalar or an array with 1-based indices.  The
result has 1-based indices and contiguous elements which is suitable for fast
linear indexing.

See also [`bcastcopy`](@ref), [`bcastdims`](@ref).

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

# Call `bcastcopy` if argument `A` can certainly not be returned.
bcastlazy(A, ::Type{T}, dims::Tuple{Vararg{Int}}) where {T} =
    bcastcopy(A, T, dims)

# Convert dimensions.
bcastlazy(A, ::Type{T}, dims::Tuple{Vararg{Integer}}) where {T} =
    bcastlazy(A, T, map(Int, dims))
bcastlazy(A, ::Type{T}, dims::Integer...) where {T} =
    bcastlazy(A, T, dims)

# Guess element type.
bcastlazy(A, dims::Tuple{Vararg{Integer}}) =
    bcastlazy(A, eltype(A), dims)
bcastlazy(A, dims::Integer...) =
    bcastlazy(A, dims)

"""

    bcastdims(size(A), size(B), ...) -> siz

yields the size `siz` of the array that would result from applying broadcasting
rules (see `broadcast`) to arguments `A`, `B`, etc.  The result is a tuple of
integers (of type `Int`).  Call [`check_dimensions`](@ref) if you want to also
make sure that the result is a list of valid dimensions.

See also [`dimensions`](@ref), [`check_dimensions`](@ref), [`bcastcopy`](@ref),
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

"""
    bcastdim(a, b) -> c

applies broadcasting rules for a single dimension (same as
`Base.Broadcasting._bcs1` but takes care of converting to `Int`), throwing an
exception if dimensions `a` and `b` are not compatible according to these
rules.

"""
bcastdim(a::Integer, b::Integer) = bcastdim(Int(a), Int(b))
bcastdim(a::Int, b::Int) =
    (a == b || b == 1 ? a : a == 1 ? b :
    throw(DimensionMismatch("arrays could not be broadcast to a common size")))
