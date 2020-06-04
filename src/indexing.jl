#
# indexing.jl --
#
# Array dimensions and indexing.
#

"""

`IndexRange` is the union of types that can define a range of indices along a
single dimension of an array, that is integers and integer-valued unit ranges.

The code purposely limited to these types, to make sure only ranges with step
equal to 1 are considered.

"""
const IndexRange = Union{<:Integer,AbstractUnitRange{<:Integer}}

# Types for the sign of an offset.
const PlusMinus = Union{typeof(+),typeof(-)}

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

See also [`has_standard_indexing`](@ref), [`same_dimensions`](@ref).

"""
dimensions(dims::Tuple{}) = dims
dimensions(dim::Integer) = (to_int(dim),)
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
    error("arrays have non-standard indices")

"""

```julia
same_dimensions(A, B...) -> size(A)
```

checks whether arrays `A`, `B`, etc., have the same dimensions and standard
indexing and returns these dimensions.  If dimensions are not all identical, a
`DimensionMismatch` exception is thrown.  If arrays have non-standard indexing
(that is indices not starting at index one), an `ArgumentError` exception is
thrown.

See also [`dimensions`](@ref), [`has_standard_indexing`](@ref).

"""
same_dimensions(A::AbstractArray) = dimensions(A)
@inline function same_dimensions(A::AbstractArray, B::AbstractArray...)
    dims = dimensions(A)
    all_match(dims, dimensions, B...) || throw_not_same_dimensions()
    return dims
end

@noinline throw_not_same_dimensions() =
    throw(DimensionMismatch("arrays must have same dimensions"))

"""

```julia
check_dimensions(dims)
```

checks whether `dims` is a list of valid dimensions.  An error is thrown if not
all dimensions are nonnegative.

See also [`dimensions`](@ref).

"""
check_dimensions(::Tuple{}) = true
check_dimensions(dims::Tuple{Vararg{Integer}}) =
    allof(nonnegative, dims...) || error("invalid array dimension(s)")
check_dimensions(dim::Integer) =
    dim ≥ 0 || error("invalid array dimension")

@deprecate checkdimensions check_dimensions

"""

```julia
nonnegative(x) -> bool
```

yields whether `x` is greater or equal zero.

"""
nonnegative(x) = x ≥ zero(x)

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
    CartesianIndices(map((i,j) -> i:j, start.I, stop.I))
end
cartesian_indices(dims::Tuple{Vararg{Integer}}) =
    CartesianIndices(map(dim -> Base.OneTo(dim), dims))
cartesian_indices(rngs::Tuple{Vararg{AbstractUnitRange{<:Integer}}}) =
    CartesianIndices(rngs)

# The following, would yield an `AbstractUnitRange` if a single argument is
# provided that is an integer or a range (not a tuple).  This lead to
# ambiguities so it is has been disabled.  The rules are that (i) an array
# argument (including a range) yields the indices for this array, (ii) a tuple
# of dimensions, a tuple of unit ranges or a pair of `CartesianIndices` yields
# the indices of the Cartesian region defined by these arguments.
#
#cartesian_indices(dim::Int) = Base.OneTo(dim)
#cartesian_indices(dim::Integer) = Base.OneTo(to_int(dim))
#cartesian_indices(rng::AbstractUnitRange{Int}) = rng
#cartesian_indices(rng::AbstractUnitRange{<:Integer}) = convert(UnitRange{Int}, rng)

"""

```julia
axis_limits(I) = (i0,i1)
```

yields the limits `i0` and `i1` of index range `I` as a 2-tuple of `Int`'s and
such that `i0:i1` represents the same indices as `I` (although not in the same
order if `step(I) < 0`).  If `step(I)` is not equal to ±1, an `ArgumentError`
exception is thrown.

"""
axis_limits(I::AbstractUnitRange{<:Integer}) =
    (to_int(first(I)), to_int(last(I)))
axis_limits(I::AbstractRange{<:Integer}) =
    ((i0, i1, s) = (to_int(first(I)), to_int(last(I)), step(I));
     (s == +1 ? (i0,i1) :
      s == -1 ? (i1,i0) : throw_invalid_range_step()))

@noinline throw_invalid_range_step() =
    throw(ArgumentError("expecting a range with a step equal to ±1"))

"""

Assuming `A` and `B` are arrays with `N` dimensions:

```julia
common_indices(A, B) -> inds
```

yields the set of all the indices that are valid for both `A` and `B`.  The
result is similar to `axes(A)` or `axes(B)`, that is an `N`-tuple of integer
valued unit ranges.

An offset `k` with a sign may be specified:

```julia
common_indices(A, B, +, k)
common_indices(A, B, -, k)
```

to obtain the set of all indices `i` such that `A[i]` and, respectively,
`B[i+k]` or `B[i-k]` are valid.  Offset `k` can be a tuple of integers or a
Cartesian index.

Arguments `A` and `B` may be both tuples of indices or index ranges or both
scalar or index range.  This is used in the following example, where we want to
do `A[i] = B[i]*C[i + k]` given the offset `k` and for all valid indices `i`:

```julia
I = common_indices(same_axes(A, B), axes(C), +, k)
@inbounds @simd for i in CartesianIndices(I)
   A[i] = B[i]*C[i + k]
end
```

Remarks: `same_axes(A,B)` is called to get the axes of `A` and `B` while
asserting that they are the same, as a result no bound checking should be
necessary and the loop can be optimzed for vectorization.

""" common_indices

@inline function common_indices(A::AbstractArray{<:Any,N},
                                B::AbstractArray{<:Any,N}) where {N}
    common_indices(axes(A), axes(B))
end
@inline function common_indices(A::AbstractArray{<:Any,N},
                                B::AbstractArray{<:Any,N},
                                op::PlusMinus,
                                K::CartesianIndex{N}) where {N}
    common_indices(axes(A), axes(B), op, K)
end
@inline function common_indices(A::AbstractArray{<:Any,N},
                                B::AbstractArray{<:Any,N},
                                op::PlusMinus,
                                K::NTuple{N,Integer}) where {N}
    common_indices(axes(A), axes(B), op, K)
end

@inline function common_indices(I::NTuple{N,IndexRange},
                                J::NTuple{N,IndexRange}) where {N}
    map(common_indices, I, J)
end

@inline function common_indices(I::NTuple{N,IndexRange},
                                J::NTuple{N,IndexRange},
                                op::PlusMinus,
                                K::CartesianIndex{N}) where {N}
    common_indices(I, J, op, K.I)
end

@inline function common_indices(I::NTuple{N,IndexRange},
                                J::NTuple{N,IndexRange},
                                ::typeof(+),
                                K::NTuple{N,Integer}) where {N}
    map(_common_indices_plus, I, J, K)
end

@inline function common_indices(I::NTuple{N,IndexRange},
                                J::NTuple{N,IndexRange},
                                ::typeof(-),
                                K::NTuple{N,Integer}) where {N}
    map(_common_indices_minus, I, J, K)
end

@inline common_indices(i::IndexRange, j::IndexRange) =
    max(to_int(first(i)), to_int(first(j))):min(to_int(last(i)), to_int(last(j)))

@inline common_indices(i::IndexRange, j::IndexRange, ::typeof(+), k::Integer) =
    _common_indices_plus(i, j, k)

@inline common_indices(i::IndexRange, j::IndexRange, ::typeof(-), k::Integer) =
    _common_indices_minus(i, j, k)

@inline _common_indices_plus(i::IndexRange, j::IndexRange, k::Integer) =
    max(to_int(first(i)), to_int(first(j) - k)):min(to_int(last(i)), to_int(last(j) - k))

@inline _common_indices_minus(i::IndexRange, j::IndexRange, k::Integer) =
    max(to_int(first(i)), to_int(first(j) + k)):min(to_int(last(i)), to_int(last(j) + k))

"""

```julia
split_interval(I, J, +/-, k) -> Ia, Ib, Ic
```

given unit ranges `I` and `J` and offset `±k`, yields 3 unit ranges, such that
`Ia ∪ Ib ∪ Ic = I` and:

- `∀ i ∈ Ia`, `i ± k < first(J)`;

- `∀ i ∈ Ib`, `i ± k ∈ J`;

- `∀ i ∈ Ic`, `i ± k > last(J)`.

unit ranges may be replaced by their first and last values:

```julia
split_interval(first(I), last(I), first(J), last(J), +/-, k)
```

yields the same result as above.

""" split_interval

@inline function split_interval(I::AbstractUnitRange{<:Integer},
                                J::AbstractUnitRange{<:Integer},
                                pm::PlusMinus, k::Integer)

    split_interval(first(I), last(I), first(J), last(J), pm, k)
end

@inline function split_interval(imin::Integer, imax::Integer,
                                jmin::Integer, jmax::Integer,
                                pm::PlusMinus, k::Integer)
    split_interval(to_int(imin), to_int(imax),
                   to_int(jmin), to_int(jmax), pm, to_int(k))
end

@inline function split_interval(imin::Int, imax::Int,
                                jmin::Int, jmax::Int, ::typeof(-), k::Int)
    split_interval(imin, imax, jmin, jmax, +, -k)
end

@inline function split_interval(imin::Int, imax::Int,
                                jmin::Int, jmax::Int, ::typeof(+), k::Int)
    #
    # Intervals for different conditions: (a) below lower bound, (b) within
    # bounds or (c) behyond upper bound.
    #
    # (a)  i ∈ imin:imax  and  j = i + k < jmin
    #      ⟺ i ∈ imin : min(imax, jmin - k - 1)
    #
    Ia = imin:min(imax,jmin-k-1)
    #
    # (b)  i ∈ imin:imax  and  jmin ≤ j = i + k ≤ jmax
    #      ⟺ i ∈ max(imin, jmin - k) ≤ i ≤ min(imax, jmax - k)
    #
    Ib = max(imin,jmin-k):min(imax,jmax-k)
    #
    # (c)  i ∈ imin:imax  and  j = i + k > jmax
    #      ⟺ i ∈ max(imin, jmax - k + 1) : imax
    #
    Ic = max(imin,jmax-k+1):imax
    return Ia, Ib, Ic
end

"""

```julia
same_axes(A, B...) -> axes(A)
```

checks whether arrays `A`, `B`, etc., have the same axes and returns them.
If axes are not all identical, a `DimensionMismatch` exception is thrown.

"""
same_axes(A::AbstractArray) = axes(A)
@inline function same_axes(A::AbstractArray, B::AbstractArray...)
    inds = axes(A)
    all_match(inds, axes, B...) || throw_not_same_axes()
    return inds
end

@noinline throw_not_same_axes() =
    throw(DimensionMismatch("arrays must have same indices"))

"""

```julia
all_indices(A...)
```

yields an iterable object for visiting each index of array(s) `A` in an
efficient manner. For array types that have opted into fast linear indexing
(like `Array`), this is simply the range `1:length(A)`. For other array types,
return a specialized Cartesian range to efficiently index into the array with
indices specified for every dimension.

If more than one `AbstractArray` argument are supplied, `all_indices` will
create an iterable object that is fast for all arguments (a `UnitRange` if all
inputs have fast linear indexing, a `CartesianIndices` otherwise).  A
`DimensionMismatch` exception is thrown if the arrays have different axes so
that it is always safe to use `@inbounds` in front of a loop like:

```julia
for i in all_indices(A, B, C, D)
   A[i] = B[i]*C[i] + D[i]
end
```

when `A`, `B` etc. are all (abstract) arrays.

This method is similar to `eachindex` except that a `DimensionMismatch`
exception is thrown if arrays have different axes.  For linearly indexed
arrays, `eachindex` only checks that they have the same linear index range
(that is the same number of elements, not the same shape).

""" all_indices

@inline all_indices(A::AbstractArray) =
    eachindex(IndexStyle(A), A)
@inline all_indices(A::AbstractArray, B::AbstractArray) =
    all_indices(IndexStyle(A, B), A, B)
@inline all_indices(A::AbstractArray, B::AbstractArray...) =
    all_indices(IndexStyle(A, B...), A, B...)

@inline function all_indices(::IndexLinear, A::AbstractArray,
                             B::AbstractArray...)
    all_match(axes(A), axes, B...) ||
        throw_indices_mismatch(IndexLinear(), A, B...)
    return eachindex(IndexLinear(), A)
end

@inline function all_indices(::IndexCartesian, A::AbstractArray,
                             B::AbstractArray...)
    inds = axes(A)
    all_match(inds, axes, B...) ||
        throw_indices_mismatch(IndexCartesian(), A, B...)
    # The following is the same as `eachindex(IndexCartesian(),A)` which yields
    # `CartesianIndices(axes(A))`.
    return CartesianIndices(inds)
end

@noinline function throw_indices_mismatch(::IndexLinear,
                                          A::AbstractArray...)
    io = IOBuffer()
    write(io, "all arguments of `all_indices` must have the same shape, got ")
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
    write(io, "all arguments of `all_indices` must have the same axes, got ")
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