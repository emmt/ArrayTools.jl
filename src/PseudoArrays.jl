#
# PseudoArrays.jl -
#
# Types and methods to facilitate the definition of custom array-like types.
#

module PseudoArrays

export
    CartesianArray,
    PseudoArray,
    LinearArray

using Base: axes1, elsize, tail, OneTo, throw_boundserror, @propagate_inbounds
import Base: getindex, setindex!, checkbounds

"""

Abstract type `PseudoArray{T,N,S}` is to be derived by types that want to
provide an array-like interface. Parameter `T` is the element type, parameter
`N` is the number of dimensions and parameter `S` is the index style:
`IndexCartesian` or `IndexLinear`.

!!! note
    The indexing style must be part of the signature because it must be
    possible to call `IndexStyle()` on the data type not the instance. Another
    possibility would have been to have the type of the embedded array be part
    of the signature but this is more restrictive.

Alias `LinearArray{T,N}` is an abstract type that can be derived by types that
want to provide an array-like interface with array values stored in an array
whose index style is linear.

Usage can be as simple as:

    struct CustomArray{T,N,...} <: LinearArray{T,N}
        arr::Array{T,N} # can be any array type with linear index style
        ...             # anything else
    end

    @inline Base.parent(A::CustomArray) = A.arr

As a result, instances of `CustomArray{T,N}` will be seen as instances of
`AbstractArray{T,N}` and behave as if they implement linear indexing. Apart
from the needs to extend the `Base.parent` method, the interface to
`LinearArray{T,N}` should provide any necessary methods for indexation, getting
the dimensions, the element type, *etc.* for the derived custom type. You may
however override these definitions by more optimized or more suitable methods
specialized for your custom array-like type.

Similarly, alias `CartesianArray{T,N}` is an abstract type that can be derived
by types that want to provide an array-like interface with array values stored
in an array whose index style is Cartesian. For such array-like object, index
checking requires an efficient implementation of the `Base.axes()` method which
you may have to specialize. The default implementation is:

    @inline Base.axes(A::PseudoArray) = axes(parent(A))

"""
abstract type PseudoArray{T,N,S<:IndexStyle} <: AbstractArray{T,N} end
const CartesianArray{T,N} = PseudoArray{T,N,IndexCartesian}
const LinearArray{T,N} = PseudoArray{T,N,IndexLinear}

@doc @doc(PseudoArray) CartesianArray
@doc @doc(PseudoArray) LinearArray

# Make PseudoArray instances behave like arrays (indexing is considered later).
#Base.eltype(::PseudoArray{T,N}) where {T,N} = T # FIXME: not needed
#Base.ndims(::PseudoArray{T,N}) where {T,N} = N # FIXME: not needed
@inline Base.length(A::PseudoArray) = length(parent(A))
@inline Base.size(A::PseudoArray) = size(parent(A))
Base.size(A::PseudoArray, d) = size(parent(A), d)
@inline Base.axes(A::PseudoArray) = axes(parent(A))
Base.axes(A::PseudoArray, d) = axes(parent(A), d)
@inline Base.axes1(A::PseudoArray) = axes1(parent(A))
@inline Base.IndexStyle(::Type{<:PseudoArray{T,N,S}}) where {T,N,S} = S()
Base.parent(A::T) where {T<:PseudoArray} =
    error(string("method parent() must be extended for instances of ", T))
Base.elsize(::Type{<:PseudoArray{T,N}}) where {T,N} = elsize(Array{T,N})
Base.sizeof(A::PseudoArray) = sizeof(parent(A))
Base.pairs(S::IndexCartesian, A::PseudoArray) = pairs(S, parent(A))
Base.pairs(S::IndexLinear, A::PseudoArray) = pairs(S, parent(A))

# Make LinearArray instances efficient iterators.
@inline Base.iterate(A::LinearArray, i=1) =
    ((i % UInt) - 1 < length(A) ? (@inbounds A[i], i + 1) : nothing)

@inline function getindex(A::LinearArray{T,N}, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds r = parent(A)[i]
    r
end

@inline function getindex(A::CartesianArray{T,N},
                          I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds r = parent(A)[I...]
    r
end

@inline function setindex!(A::LinearArray{T,N}, val, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i] = val
    A
end

@inline function setindex!(A::CartesianArray{T,N}, val,
                           I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds parent(A)[I...] = val
    A
end

end # module
