"""

```julia
rubberindex(n)
```

yields a rubber index of lenght `n`.  That is a `n`-tuple of colons `:`.

When `n` is known at compile time, it is faster to call:

```julia
rubberindex(Val(n))
```

This method is suitable to extract sub-arrays of build views when some kind of
rubber index is needed.  For instance:

```julia
slice(A::AbstractArray{T,N}, i::Integer) where {T,N} =
    A[rubberindex(Val(N-1))..., i]
```

defines a function that returns the `i`-th slice of `A` assuming index `i`
refers the last index of `A`.

See also: [`…`](@ref), [`RubberIndex`](@ref).

"""
rubberindex(n::Integer) =
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
     _rubberindex(n))

function _rubberindex(n::Integer)
    n ≥ 0 || throw(ArgumentError(string("number of dimensions should be ≥ 0, got ", n)))
    return ([Colon() for i in 1:n]...,)
end

rubberindex(::Val{ 0}) = ()
rubberindex(::Val{ 1}) = (:,)
rubberindex(::Val{ 2}) = (:,:,)
rubberindex(::Val{ 3}) = (:,:,:,)
rubberindex(::Val{ 4}) = (:,:,:,:,)
rubberindex(::Val{ 5}) = (:,:,:,:,:,)
rubberindex(::Val{ 6}) = (:,:,:,:,:,:,)
rubberindex(::Val{ 7}) = (:,:,:,:,:,:,:,)
rubberindex(::Val{ 8}) = (:,:,:,:,:,:,:,:,)
rubberindex(::Val{ 9}) = (:,:,:,:,:,:,:,:,:,)
rubberindex(::Val{10}) = (:,:,:,:,:,:,:,:,:,:,)
rubberindex(v::Val{N}) where {N} = ntuple(x -> :, v)

struct RubberIndex; end

"""

`RubberIndex` is the singleron type that represents any number of indices.  The
constant `…` is defined as `RubberIndex()` and can be used in array indexation
to left and/or right justify the other indices.  For instance, assuming `A` is
a `3×4×5×6` array, then all the following equalities hold:

```julia
A[…]           == A[:,:,:,:]
A[…,3]         == A[:,:,:,3]
A[2,…]         == A[2,:,:,:]
A[…,2:4,5]     == A[:,:,2:4,5]
A[2:3,…,1,2:4] == A[2:3,:,1,2:4]
```

As you can see, the advantage of the rubber index `…` is that it automatically
expands as the number of colons needed to have the correct number of indices.
The expressions are also more readable.

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

!!! warn
    There is two known limitation:
    1. The `end` reserved word can only be used in intervals specified *before*
       the rubber index but not *after*.  This limitation is due to the Julia
       parser cannot be avoided.
    2. At most 9 indices can be specified before the rubber index.  This
       can be extended by editing the source code.

See also: [`rubberindex`](@ref).

"""
const … = RubberIndex()
@doc @doc(…) RubberIndex

"""

`Index` is the union of types that are eligible as a single index, that is
integers and integer valued ranges.

"""
const Index = Union{Integer,AbstractRange{<:Integer},Colon,CartesianIndex}

"""

```julia
numberofindices(inds...)
```

yields the total number of indices specified by `inds...`.  Integers, colons
and integer valued ranges count as one each, Cartesian indices count as their
dimensionality.

"""
numberofindices() = 0
numberofindices(::Colon) = 1
numberofindices(::Integer) = 1
numberofindices(::CartesianIndex{N}) where {N} = N
numberofindices(::AbstractRange{<:Integer}) = 1
numberofindices(i::Index, inds::Index...) = numberofindices(i) + numberofindices(inds...)

# Generate a tuple with as many colons as the number of dimensions `N` minus
# the number of indices specified by `inds...`.
_colons(N::Int, inds::Index...) = begin
    (n = N - numberofindices(inds...)) ≥ 0 ||
         throw(DimensionMismatch("the number of specified indices exceed the number of dimensions"))
    rubberindex(n)
end

# A[…]
getindex(A::AbstractArray, ::RubberIndex) = copy(A)
setindex!(A::AbstractArray{T,N}, val, ::RubberIndex) where {T,N} =
    A[rubberindex(Val(N))...] = val
dotview(A::AbstractArray{T,N}, ::RubberIndex) where {T,N} =
    dotview(A, rubberindex(Val(N))...)

# A[…, inds...]
function getindex(A::AbstractArray{T,N},
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[_colons(N, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[_colons(N, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, _colons(N, inds...)..., inds...)
end

# A[i, …]
function getindex(A::AbstractArray{T,N},
                  i::Index,
                  ::RubberIndex) where {T,N}
    A[i, _colons(N, i)...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i::Index,
                   ::RubberIndex) where {T,N}
    A[i, _colons(N, i)...] = val
end
function dotview(A::AbstractArray{T,N},
                 i::Index,
                 ::RubberIndex) where {T,N}
    dotview(A, i, _colons(N, i)...)
end

# A[i, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i, _colons(N, i, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i, _colons(N, i, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i, _colons(N, i, inds...)..., inds...)
end

# A[i1, i2, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, _colons(N, i1, i2, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, _colons(N, i1, i2, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, _colons(N, i1, i2, inds...)..., inds...)
end

# A[i1, i2, i3, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, _colons(N, i1, i2, i3, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, _colons(N, i1, i2, i3, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, _colons(N, i1, i2, i3, inds...)..., inds...)
end

# A[i1, i2, i3, i4, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index, i4::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, i4, _colons(N, i1, i2, i3, i4, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index, i4::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, i4, _colons(N, i1, i2, i3, i4, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index, i4::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, i4, _colons(N, i1, i2, i3, i4, inds...)..., inds...)
end

# A[i1, i2, i3, i4, i5, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5,
      _colons(N, i1, i2, i3, i4, i5, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5,
      _colons(N, i1, i2, i3, i4, i5, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, i4, i5,
            _colons(N, i1, i2, i3, i4, i5, inds...)..., inds...)
end

# A[i1, i2, i3, i4, i5, i6, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                  i6::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6,
      _colons(N, i1, i2, i3, i4, i5, i6, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                   i6::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6,
      _colons(N, i1, i2, i3, i4, i5, i6, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                 i6::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, i4, i5, i6,
            _colons(N, i1, i2, i3, i4, i5, i6, inds...)..., inds...)
end

# A[i1, i2, i3, i4, i5, i6, i7, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                  i6::Index, i7::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6, i7,
      _colons(N, i1, i2, i3, i4, i5, i6, i7, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                   i6::Index, i7::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6, i7,
      _colons(N, i1, i2, i3, i4, i5, i6, i7, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                 i6::Index, i7::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, i4, i5, i6, i7,
            _colons(N, i1, i2, i3, i4, i5, i6, i7, inds...)..., inds...)
end

# A[i1, i2, i3, i4, i5, i6, i7, i8, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                  i6::Index, i7::Index, i8::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6, i7, i8,
      _colons(N, i1, i2, i3, i4, i5, i6, i7, i8, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                   i6::Index, i7::Index, i8::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6, i7, i8,
      _colons(N, i1, i2, i3, i4, i5, i6, i7, i8, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                 i6::Index, i7::Index, i8::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, i4, i5, i6, i7, i8,
            _colons(N, i1, i2, i3, i4, i5, i6, i7, i8, inds...)..., inds...)
end

# A[i1, i2, i3, i4, i5, i6, i7, i8, i9, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                  i6::Index, i7::Index, i8::Index, i9::Index,
                  ::RubberIndex,
                  inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6, i7, i8, i9,
      _colons(N, i1, i2, i3, i4, i5, i6, i7, i8, i9, inds...)..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                   i6::Index, i7::Index, i8::Index, i9::Index,
                   ::RubberIndex,
                   inds::Index...) where {T,N}
    A[i1, i2, i3, i4, i5, i6, i7, i8, i9,
      _colons(N, i1, i2, i3, i4, i5, i6, i7, i8, i9, inds...)..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Index, i2::Index, i3::Index, i4::Index, i5::Index,
                 i6::Index, i7::Index, i8::Index, i9::Index,
                 ::RubberIndex,
                 inds::Index...) where {T,N}
    dotview(A, i1, i2, i3, i4, i5, i6, i7, i8, i9,
            _colons(N, i1, i2, i3, i4, i5, i6, i7, i8, i9, inds...)..., inds...)
end
