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

"""
    PseudoArray{T,N,S}

is the type abstract to be derived by types that wrap a *parent* array-like
object and extend the abstract array API to access this object. Parameter `T`
is the element type, parameter `N` is the number of dimensions and parameter
`S` is the index style: `IndexCartesian` or `IndexLinear`.

Usage can be as simple as:

    struct CustomArray{T,N,...} <: PseudoArray{T,N,IndexLinear}
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

The following base methods are extended for pseudo-arrays `A` of type
`PseudoArray{T,N,S}`:

- `Base.length(A)`, `Base.size(A[,d])`, `Base.axes(A[,d])`, and `Base.axes1(A)`.

- `Base.getindex(A,i)` and `Base.setindex!(A,i)` for `i` of type `Int` if
  `S` is `IndexLinear`, of type `Vararg{Int,N}` otherwise.

- `Base.iterate(A,state)` if `S` is `IndexLinear`; otherwise the default
  definition of `Base.iterate` is used.


# See also

[`CartesianArray`](@ref) and [`LinearArray`](@ref).

"""
abstract type PseudoArray{T,N,S<:IndexStyle} <: AbstractArray{T,N} end

"""
    const CartesianArray{T,N} = PseudoArray{T,N,IndexCartesian}

is an alias to pseudo-arrays with Cartesian indexing.

# See also

[`PseudoArray`](@ref) and [`LinearArray`](@ref).

"""
const CartesianArray{T,N} = PseudoArray{T,N,IndexCartesian}

"""
    const LinearArray{T,N} = PseudoArray{T,N,IndexLinear}

is an alias to pseudo-arrays with linear indexing.

# See also

[`PseudoArray`](@ref) and [`CartesianArray`](@ref).

"""
const LinearArray{T,N} = PseudoArray{T,N,IndexLinear}

# Make PseudoArray instances behave like arrays (indexing is considered later).
@inline Base.length(A::PseudoArray) = length(parent(A))
@inline Base.size(A::PseudoArray) = size(parent(A))
Base.size(A::PseudoArray, d) = size(parent(A), d)
@inline Base.axes(A::PseudoArray) = axes(parent(A))
Base.axes(A::PseudoArray, d) = axes(parent(A), d)
@inline Base.axes1(A::PseudoArray) = Base.axes1(parent(A))
@inline Base.IndexStyle(::Type{<:PseudoArray{T,N,S}}) where {T,N,S} = S()
Base.parent(A::T) where {T<:PseudoArray} =
    error("method `parent(A::T)` must be extended for instances of `$T`")
Base.elsize(::Type{<:PseudoArray{T,N}}) where {T,N} = Base.elsize(Array{T,N})
Base.sizeof(A::PseudoArray) = sizeof(parent(A))
Base.pairs(S::IndexCartesian, A::PseudoArray) = pairs(S, parent(A))
Base.pairs(S::IndexLinear, A::PseudoArray) = pairs(S, parent(A))

# Make LinearArray instances efficient iterators.
@inline Base.iterate(A::LinearArray, i=1) =
    ((i % UInt) - 1 < length(A) ? (@inbounds A[i], i + 1) : nothing)

@inline function Base.getindex(A::LinearArray{T,N}, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    return @inbounds parent(A)[i]
end

@inline function Base.getindex(A::CartesianArray{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    return @inbounds parent(A)[I...]
end

@inline function Base.setindex!(A::LinearArray{T,N}, val, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i] = val
    return A
end

@inline function Base.setindex!(A::CartesianArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds parent(A)[I...] = val
    return A
end

end # module
