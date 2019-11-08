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
to left or right justify the other indices.  For instance, assuming `A` is a
`3×4×5×6` array, then all the following equalities hold:

```julia
A[…]           == A[:,:,:,:]
A[…,3]         == A[:,:,:,3]
A[2,…]         == A[2,:,:,:]
A[…,2:4,5]     == A[:,:,2:4,5]
A[2:3,…,1,2:4] == A[2:3,:,1,2:4]
```

As you can see the advantage of the rubber index `…` is that it automatically
expands as the number of colons needed to have the correct number of indices.
The expressions are also more readable.

The rubber index may also be used for setting values.  For instance:

```julia
A[…] .= 1       # to fill A with ones
A[…,3] = A[…,2] # to copy A[:,:,:,2] in A[:,:,:,3]
A[2,…] = A[3,…] # to copy A[3,:,:,:] in A[2,:,:,:]
A[…,2:4,5] .= 7 # to set all elements in A[:,:,2:4,5] to 7
```

Leading/trailing indices may be specified as Cartesian indices (of type
`CartesianIndex`).

See also: [`rubberindex`](@ref).

"""
const … = RubberIndex()
@doc @doc(…) RubberIndex

"""

`Indices` is the union of types that are eligible as a single index, that is
integers and integer valued ranges.

"""
const Indices = Union{Integer,AbstractRange{<:Integer}}

#numberofindices(::Integer) = 1
#numberofindices(::CartesianIndex{N}) where {N} = N
#numberofindices(::AbstractRange{<:Integer}) = 1
#numberofindices(inds::Indices...) = length(args)

# A[…]
getindex(A::AbstractArray, ::RubberIndex) = copy(A)
setindex!(A::AbstractArray{T,N}, val, ::RubberIndex) where {T,N} =
    A[rubberindex(Val(N))] = val
dotview(A::AbstractArray{T,N}, ::RubberIndex) where {T,N} =
    dotview(A, rubberindex(Val(N))...)

# A[…, i]
function getindex(A::AbstractArray{T,N},
                  ::RubberIndex,
                  i::Indices) where {T,N}
    A[rubberindex(N - 1)..., i]
end
function setindex!(A::AbstractArray{T,N}, val,
                   ::RubberIndex,
                   i::Indices) where {T,N}
    A[rubberindex(N - 1)..., i] = val
end
function dotview(A::AbstractArray{T,N},
                 ::RubberIndex,
                 i::Indices) where {T,N}
    dotview(A, rubberindex(N - 1)..., i)
end

# A[…, inds...]
function getindex(A::AbstractArray{T,N},
                  ::RubberIndex,
                  inds::Indices...) where {T,N}
    A[rubberindex(N - length(inds))..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   ::RubberIndex,
                   inds::Indices...) where {T,N}
    A[rubberindex(N - length(inds))..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 ::RubberIndex,
                 inds::Indices...) where {T,N}
    dotview(A, rubberindex(N - length(inds))..., inds...)
end

# A[i, …]
function getindex(A::AbstractArray{T,N},
                  i::Indices,
                  ::RubberIndex) where {T,N}
    A[i, rubberindex(N - 1)...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i::Indices,
                   ::RubberIndex) where {T,N}
    A[i, rubberindex(N - 1)...] = val
end
function dotview(A::AbstractArray{T,N},
                 i::Indices,
                 ::RubberIndex) where {T,N}
    dotview(A, i, rubberindex(N - 1)...)
end

# A[i, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i::Indices,
                  ::RubberIndex,
                  inds::Indices...) where {T,N}
    A[i, rubberindex(N - 1 - length(inds))..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i::Indices,
                   ::RubberIndex,
                   inds::Indices...) where {T,N}
    A[i, rubberindex(N - 1 - length(inds))..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i::Indices,
                 ::RubberIndex,
                 inds::Indices...) where {T,N}
    dotview(A, i, rubberindex(N - 1 - length(inds))..., inds...)
end

# A[i1, i2, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Indices,
                  i2::Indices,
                  ::RubberIndex,
                  inds::Indices...) where {T,N}
    A[i1, i2, rubberindex(N - 2 - length(inds))..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Indices,
                   i2::Indices,
                   ::RubberIndex,
                   inds::Indices...) where {T,N}
    A[i1, i2, rubberindex(N - 2 - length(inds))..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Indices,
                 i2::Indices,
                 ::RubberIndex,
                 inds::Indices...) where {T,N}
    dotview(A, i1, i2, rubberindex(N - 2 - length(inds))..., inds...)
end

# A[i1, i2, i3, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Indices,
                  i2::Indices,
                  i3::Indices,
                  ::RubberIndex,
                  inds::Indices...) where {T,N}
    A[i1, i2, i3, rubberindex(N - 3 - length(inds))..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Indices,
                   i2::Indices,
                   i3::Indices,
                   ::RubberIndex,
                   inds::Indices...) where {T,N}
    A[i1, i2, i3, rubberindex(N - 3 - length(inds))..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Indices,
                 i2::Indices,
                 i3::Indices,
                 ::RubberIndex,
                 inds::Indices...) where {T,N}
    dotview(A, i1, i2, i3, rubberindex(N - 3 - length(inds))..., inds...)
end

# A[i1, i2, i3, i4, …, inds...]
function getindex(A::AbstractArray{T,N},
                  i1::Indices,
                  i2::Indices,
                  i3::Indices,
                  i4::Indices,
                  ::RubberIndex,
                  inds::Indices...) where {T,N}
    A[i1, i2, i3, i4, rubberindex(N - 4 - length(inds))..., inds...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   i1::Indices,
                   i2::Indices,
                   i3::Indices,
                   i4::Indices,
                   ::RubberIndex,
                   inds::Indices...) where {T,N}
    A[i1, i2, i3, i4, rubberindex(N - 4 - length(inds))..., inds...] = val
end
function dotview(A::AbstractArray{T,N},
                 i1::Indices,
                 i2::Indices,
                 i3::Indices,
                 i4::Indices,
                 ::RubberIndex,
                 inds::Indices...) where {T,N}
    dotview(A, i1, i2, i3, i4, rubberindex(N - 4 - length(inds))..., inds...)
end

# A[…, I]
function getindex(A::AbstractArray{T,N},
                  ::RubberIndex,
                  I::CartesianIndex{L}) where {T,N,L}
    A[rubberindex(N - L)..., I]
end
function setindex!(A::AbstractArray{T,N}, val,
                   ::RubberIndex,
                   I::CartesianIndex{L}) where {T,N,L}
    A[rubberindex(N - L)..., I] = val
end
function dotview(A::AbstractArray{T,N},
                 ::RubberIndex,
                 I::CartesianIndex{L}) where {T,N,L}
    dotview(A, rubberindex(N - L)..., I)
end

# A[I, …]
function getindex(A::AbstractArray{T,N},
                  I::CartesianIndex{L},
                  ::RubberIndex) where {T,N,L}
    A[I, rubberindex(N - L)...]
end
function setindex!(A::AbstractArray{T,N}, val,
                   I::CartesianIndex{L},
                   ::RubberIndex) where {T,N,L}
    A[I, rubberindex(N - L)...] = val
end
function dotview(A::AbstractArray{T,N},
                 I::CartesianIndex{L},
                 ::RubberIndex) where {T,N,L}
    dotview(A, I, rubberindex(N - L)...)
end

# A[I1, …, I2]
function getindex(A::AbstractArray{T,N},
                  I1::CartesianIndex{L1},
                  ::RubberIndex,
                  I2::CartesianIndex{L2}) where {T,N,L1,L2}
    A[I1, rubberindex(N - L1 - L2)..., I2]
end
function setindex!(A::AbstractArray{T,N}, val,
                   I1::CartesianIndex{L1},
                   ::RubberIndex,
                   I2::CartesianIndex{L2}) where {T,N,L1,L2}
    A[I1, rubberindex(N - L1 - L2)..., I2] = val
end
function dotview(A::AbstractArray{T,N},
                 I1::CartesianIndex{L1},
                 ::RubberIndex,
                 I2::CartesianIndex{L2}) where {T,N,L1,L2}
    dotview(A, I1, rubberindex(N - L1 - L2)..., I2)
end
