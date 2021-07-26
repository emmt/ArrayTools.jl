#
# rubberindex.jl -
#
# Rubber-index for Julia arrays.
#

"""
    colons(n)

yields a `n`-tuple of colons `:` (a.k.a. `Colon()`).

When `n` is known at compile time, it is faster to call:

    colons(Val(n))

This method is suitable to extract sub-arrays of build views when some kind of
rubber index is needed.  For instance:

    slice(A::AbstractArray{T,N}, i::Integer) where {T,N} =
        A[colons(Val(N-1))..., i]

defines a function that returns the `i`-th slice of `A` assuming index `i`
refers the last index of `A`.  Using the rubber-index `..`, a shorter
definition is:

    slice(A::AbstractArray, i) = A[.., i]

which is also able to deal with multiple trailing indices if `i` is a
`CartesianIndex`.

See also: `..`, [`RubberIndex`](@ref).

"""
colons(n::Integer) =
    (n ==  0 ? () :
     n ==  1 ? (:,) :
     n ==  2 ? (:,:,) :
     n ==  3 ? (:,:,:,) :
     n ==  4 ? (:,:,:,:,) :
     n ==  5 ? (:,:,:,:,:,) :
     n ==  6 ? (:,:,:,:,:,:,) :
     n ==  7 ? (:,:,:,:,:,:,:,) :
     n ==  8 ? (:,:,:,:,:,:,:,:,) :
     n ==  9 ? (:,:,:,:,:,:,:,:,:,) :
     n == 10 ? (:,:,:,:,:,:,:,:,:,:,) :
     _colons(n))

function _colons(n::Integer)
    n ≥ 0 || bad_ndims(n)
    return ([Colon() for i in 1:n]...,)
end

@noinline bad_ndims(n::Integer) =
    throw(ArgumentError(string("number of dimensions should be ≥ 0, got ", n)))

colons(v::Val{N}) where {N} = ntuple(colon, v)

# Just yields a colon whatever the argument.
colon(x) = Colon()

"""

`RubberIndex` is the singleron type that represents any number of indices.  The
constant `..` is defined as `RubberIndex()` and can be used in array indexation
to left and/or right justify the other indices.  For instance, assuming `A` is
a `3×4×5×6` array, then all the following equalities hold:

    A[..]           == A[:,:,:,:]
    A[..,3]         == A[:,:,:,3]
    A[2,..]         == A[2,:,:,:]
    A[..,2:4,5]     == A[:,:,2:4,5]
    A[2:3,..,1,2:4] == A[2:3,:,1,2:4]

As you can see, the advantage of the rubber index `..` is that it automatically
expands as the number of colons needed to have the correct number of indices.
The expressions are also more readable.

The rubber index may also be used for setting values.  For instance:

    A[..] .= 1         # to fill A with ones
    A[..,3] = A[..,2]  # to copy A[:,:,:,2] in A[:,:,:,3]
    A[..,3] .= A[..,2] # idem but faster
    A[2,..] = A[3,..]  # to copy A[3,:,:,:] in A[2,:,:,:]
    A[..,2:4,5] .= 7   # to set all elements in A[:,:,2:4,5] to 7

Leading/trailing indices may be specified as Cartesian indices (of type
`CartesianIndex`).

!!! warning
    There are two known limitations:
    1. The `end` reserved word can only be used in intervals specified *before*
       the rubber index but not *after*.  This limitation is due to the Julia
       parser cannot be avoided.
    2. At most 9 indices can be specified before the rubber index.  This
       can be extended by editing the source code.

See also: [`colons`](@ref).

""" RubberIndex

struct RubberIndex end

const .. = RubberIndex()
const … = .. # FIXME: should be deprecated

# Quickly get the tuple inside a Cartesian index.
to_tuple(I::CartesianIndex) = I.I

# Grow tuple of colons until tuple of axes is empty.
@inline growcolons(colons, inds::Tuple{}) = colons
@inline growcolons(colons, inds::Tuple) = growcolons((colons..., :), tail(inds))

# Drop as many colons as there are specified indices.
@inline dropcolons(colons, I::Tuple{}) = colons
@inline dropcolons(colons::Tuple{}, ::Tuple{}) = ()
@noinline dropcolons(colons::Tuple{}, I::Tuple) =
    throw(ArgumentError("too many indices specifed"))
@inline dropcolons(colons, I::Tuple) =
    dropcolons(tail(colons), tail(I))
@inline dropcolons(colons, I::Tuple{CartesianIndex, Vararg}) =
    dropcolons(dropcolons(colons, to_tuple(I[1])), tail(I))
@noinline dropcolons(colons, I::Tuple{RubberIndex, Vararg}) =
    throw(ArgumentError("more than one rubber index specified"))

@inline function to_indices(A, inds, I::Tuple{RubberIndex, Vararg})
    # Align the remaining indices to the tail of the `inds`.  First
    # `growcolons` is called to build a tuple of as many as colons as the
    # number of remaining axes. Second, `dropcolons` is used to drop as
    # many colons as the number of specified indices (taking into account
    # Cartesian indices).
    colons = dropcolons(growcolons((), inds), tail(I))
    to_indices(A, inds, (colons..., tail(I)...))
end

# avoid copying if indexing with .. alone, see
# https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/214
@inline Base.getindex(A::AbstractArray, ::RubberIndex) = A

# The following is needed to allow for statements like `A[..] .+= expr` to
# work properly.
dotview(A::AbstractArray{T,N}, ::RubberIndex) where {T,N} =
    dotview(A, ntuple(colon, Val{N}())...)
