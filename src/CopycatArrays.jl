#
# CopycatArrays.jl -
#
# Types and methods to facilitate the definition of custom array-like types.
#
#-------------------------------------------------------------------------------
#
# This file if part of the TAO software (https://github.com/emmt/ArrayTools)
# licensed under the MIT license.
#
# Copyright (C) 2019, Éric Thiébaut.
#

module CopycatArrays

export
    CartesianArray,
    CopycatArray,
    LinearArray

using Base: axes1, elsize, tail, OneTo, throw_boundserror, @propagate_inbounds
import Base: getindex, setindex!, checkbounds

"""

Abstract type `CopycatArray{T,N,S}` is to be derived by types that want to
provide an array-like interface.  Parameter `T` is the element type, parameter
`N` is the number of dimensions and parameter `S` is the index style:
`IndexCartesian` or `IndexLinear`.


!!! note
    The indexing style must be part of the signature because it must be
    possible to call `IndexStyle()` on the data type not the instance.  Another
    possibility would have been to have the type of the embedded array be part
    of the signature but this is more restrictive.

Alias `LinearArray{T,N}` is an abstract type that can be derived by types that
want to provide an array-like interface with array values stored in an array
whose index style is linear.

Usage can be as simple as:

```julia
struct CustomArray{T,N,...} <: LinearArray{T,N}
    arr::Array{T,N} # can be any array type with linear index style
    ...             # anything else
end

@inline Base.parent(A::CustomArray) = A.arr
```

As a result, instances of `CustomArray{T,N}` will be seen as instances of
`AbstractArray{T,N}` and behave as if they implement linear indexing.  Apart
from the needs to extend the `Base.parent` method, the interface to
`LinearArray{T,N}` should provide any necessary methods for indexation, getting
the dimensions, the element type, *etc.* for the derived custom type.  You may
however override these definitions by more optimized or more suitable methods
specialized for your custom array-like type.

Similarly, alias `CartesianArray{T,N}` is an abstract type that can be derived
by types that want to provide an array-like interface with array values stored
in an array whose index style is Cartesian.  For such array-like object, index
checking requires an efficient implementation of the `Base.axes()` method which
you may have to specialize.  The default implementation is:

```julia
@inline Base.axes(A::CopycatArray) = axes(parent(A))
```

"""
abstract type CopycatArray{T,N,S<:IndexStyle} <: AbstractArray{T,N} end
const CartesianArray{T,N} = CopycatArray{T,N,IndexCartesian}
const LinearArray{T,N} = CopycatArray{T,N,IndexLinear}

@doc @doc(CopycatArray) CartesianArray
@doc @doc(CopycatArray) LinearArray

# Make CopycatArray instances behave like arrays (indexing is considered later).
#Base.eltype(::CopycatArray{T,N}) where {T,N} = T # FIXME: not needed
#Base.ndims(::CopycatArray{T,N}) where {T,N} = N # FIXME: not needed
@inline Base.length(A::CopycatArray) = length(parent(A))
@inline Base.size(A::CopycatArray) = size(parent(A))
Base.size(A::CopycatArray, d) = size(parent(A), d)
@inline Base.axes(A::CopycatArray) = axes(parent(A))
Base.axes(A::CopycatArray, d) = axes(parent(A), d)
@inline Base.axes1(A::CopycatArray) = axes1(parent(A))
@inline Base.IndexStyle(::Type{<:CopycatArray{T,N,S}}) where {T,N,S} = S()
Base.parent(A::T) where {T<:CopycatArray} =
    error(string("method parent() must be extended for instances of ", T))
Base.elsize(::Type{<:CopycatArray{T,N}}) where {T,N} = elsize(Array{T,N})
Base.sizeof(A::CopycatArray) = sizeof(parent(A))
Base.pairs(S::IndexCartesian, A::CopycatArray) = pairs(S, parent(A))
Base.pairs(S::IndexLinear, A::CopycatArray) = pairs(S, parent(A))

# Make LinearArray instances efficient iterators.
@inline Base.iterate(A::LinearArray, i=1) =
    ((i % UInt) - 1 < length(A) ? (@inbounds A[i], i + 1) : nothing)

@inline @propagate_inbounds function getindex(A::LinearArray{T,N},
                                              i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds r = getindex(parent(A), i)
    return r
end

@inline @propagate_inbounds function getindex(A::CartesianArray{T,N},
                                              I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds r = getindex(parent(A), I...)
    return r
end

@inline @propagate_inbounds function setindex!(A::LinearArray{T,N}, x,
                                               i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds r = setindex!(parent(A), x, i)
    return r
end

@inline @propagate_inbounds function setindex!(A::CartesianArray{T,N}, x,
                                               I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds r = setindex!(parent(A), x, I...)
    return r
end

# FIXME: not needed?
@inline checkbounds(A::LinearArray{T,N}, i::Int) where {T,N} =
    1 ≤ i ≤ length(A) || throw_boundserror(A, i)

end # module
