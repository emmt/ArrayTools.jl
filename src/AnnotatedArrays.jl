"""

Module `AnnotatedArrays` implements annotated arrays which combine values stored in an array
and properties stored in a dictionary or a named tuple.

"""
module AnnotatedArrays

export
    AnnotatedArray,
    properties

using TypeUtils: @public
@public(
    AbstractAnnotatedArray,
    is_key)

using Base: @propagate_inbounds

#using ..ArrayTools.PseudoArrays

using TypeUtils

# Private singleton used to indicate an entry which is not found.
struct NotFound end

# Union of callable for do-block syntax.
const Callable = Union{Function, Type}

# Union of allowed types to store the properties of annotated arrays.
const Properties = Union{AbstractDict,NamedTuple}

# Union of types to storing the properties of annotated arrays with symbolic keys.
const SymbolicProperties = Union{AbstractDict{Symbol},NamedTuple}

"""
    ArrayTools.AnnotatedArrays.AbstractAnnotatedArray{T,N,P}

Abstract super-type of annotated arrays which store values of type `T` as in a
`N`-dimensional array with indexing style `S` and properties in an object of type `P` (a
dictionary or a named tuple).

"""
abstract type AbstractAnnotatedArray{T,N,P<:Properties} <: AbstractArray{T,N} end

struct AnnotatedArray{T,N,
                      P<:Properties,
                      A<:AbstractArray{T,N}} <: AbstractAnnotatedArray{T,N,P}
    array::A
    props::P

    # Private inner constructor. Should only be called with checked arguments.
    global _AnnotatedArray
    function _AnnotatedArray(array::A, props::P) where {T,N,
                                                        A<:AbstractArray{T,N},
                                                        P<:Union{NamedTuple,AbstractDict}}
        return new{T,N,P,A}(array, props)
    end
end

# Accessor for the array entries.
@inline Base.values(A::AnnotatedArray) = getfield(A, :array)

# Abstract array API for annotated arrays.
Base.length(A::AnnotatedArray) = length(values(A))
Base.size(A::AnnotatedArray) = size(values(A))
Base.axes(A::AnnotatedArray) = axes(values(A))
Base.IndexStyle(::Type{<:AnnotatedArray{T,N,P,A}}) where {T,N,P,A} = IndexStyle(A)

"""
    A = AnnotatedArray{T}(undef, dims, props...)
    A = AnnotatedArray{T}(undef, dims...; kwds...)

Return an annotated array storing its values in a new array with undefined elements of type
`T` and size `dims` and with properties specified by `props...` (a dictionary, a named
tuple, or any number of `key=>val` pairs), or by the keywords `kwds...`.

If properties are specified as a named tuple, the properties are immutable; otherwise,
properties are mutable. The properties of the annotated array are indexed by symbols if they
are initially specified by keywords, by a named tuple, or by a dictionary or pairs with
symbolic keys. Otherwise, the key type of the properties corresponds to that of the
dictionary or pairs.

If initial properties are specified by a dictionary, this dictionary is shared by the
annotated array.

For symbolic keys, accessing or mutating a property of the annotated array `A` is done by
`A.key` and `A.key = val`.

Expressions `A[i]` and `A[i] = val` can be used to access or mutate an entry (if `i` is an
*index*) or a property (if `i` is a *key*) of the annotated array `A`. However, to avoid
confusing property keys and array indices, `i` is only considered as being a property key if
it is a *string* or a *symbol*. The fucntion [`AnnotatedArray.is_key`](@ref
ArrayTools.AnnotatedArray.is_key) may be specialized to overcome this rule if necessary.
Other base methods such as `get`, `get!`, and `hashey` also rely on `AnnotatedArray.is_key`
to disentangle between property keys and array indices.

Call [`properties(A)`](@ref ArrayTools.AnnotatedArray.properties) or `values(A)` to retrieve
the respective objects storing the properties or the array entries of the an annotated array
`A`. For example, the number of properties is given by `length(properties(A))` while the
number of entries is given by `length(values(A))` or by `length(A)`.

"""
AnnotatedArray{T,N}(::UndefInitializer, dims::Integer...; kwds...) where {T,N} =
    AnnotatedArray{T,N}(undef, dims; kwds...)
AnnotatedArray{T}(::UndefInitializer, dims::Integer...; kwds...) where {T} =
    AnnotatedArray{T}(init, dims; kwds...)
AnnotatedArray{T,N}(::UndefInitializer, dims::NTuple{N,Integer}, args...; kwds...) where {T,N} =
    AnnotatedArray{T}(undef, dims, args...; kwds...)
AnnotatedArray{T}(::UndefInitializer, dims::Tuple{Vararg{Integer}}, args...; kwds...) where {T} =
    AnnotatedArray(Array{T}(undef, dims), args...; kwds...)

"""
    A = AnnotatedArray(arr::AbstractArray, props...)
    A = AnnotatedArray(arr::AbstractArray; kwds...)

Return an annotated array storing its values in array `arr` and with initial properties
specified by `props...` or by `kwds...`.

"""
AnnotatedArray(array::AbstractArray, props::NamedTuple) = _AnnotatedArray(array, props)

function AnnotatedArray(array::AbstractArray, props::AbstractDict{K,V}) where {K,V}
    K === Any && error("key type must be more specialized than `Any`")
    K <: Integer && error("key type must not be an integer")
    K <: CartesianIndex && error("key type must not be a Cartesian index")
    K <: AbstractRange{<:Integer} && error("key type must not be an integer range")
    K <: Colon && error("key type must not be a colon")
    return _AnnotatedArray(array, props)
end

# Get rid of matching type-parameter when initial array is specified. The pairs of
# definitions are needed to disentangle ambiguities.
AnnotatedArray{T,N}(array::AbstractArray{T,N}, args...; kwds...) where {T,N} =
    AnnotatedArray(array, args...; kwds...)
#AnnotatedArray{T,N}(array::AbstractArray{T,N}; kwds...) where {T,N} =
#    AnnotatedArray(array; kwds...)
AnnotatedArray{T}(array::AbstractArray{T}, args...; kwds...) where {T} =
    AnnotatedArray(array, args...; kwds...)
#AnnotatedArray{T}(array::AbstractArray{T}; kwds...) where {T} =
#    AnnotatedArray(array; kwds...)
AnnotatedArray{T}(array::AbstractArray, args...; kwds...) where {T} =
    AnnotatedArray(convert_eltype(T, array)::AbstractArray{T}, args...; kwds...)

# Constructor for properties provided by keywords.
AnnotatedArray(array::AbstractArray; kwds...) =
    AnnotatedArray(array, Dict{Symbol,Any}(kwds...))

# Constructors with initial array and initial properties specified as key-value pairs. When
# initial properties are specified as key-value pairs, we want to have a dictionary whose
# key type is specialized if possible (for efficiency) but avoid having value type
# specialized (for flexibility). If one wants specific properties key and value types, it is
# always possible to explicitly specify a dictionary.
function AnnotatedArray(array::AbstractArray, p::Pair, ps::Pair...)
    K = infer_keytype(keytype(p), ps...)
    return AnnotatedArray(array, Dict{K,Any}(p, ps...))
end

infer_keytype(::Type{T}) where {T} = T
@inline infer_keytype(::Type{T}, p::Pair, ps::Pair...) where {T} =
    infer_keytype(promote_keytype(T, keytype(p)), ps...)

promote_keytype(::Type{<:AbstractString}, ::Type{<:AbstractString}) = String
@inline function promote_keytype(::Type{S}, ::Type{T}) where {S,T}
    S >: T && return S
    S <: T && return T
    return promote_keytype(supertype(S), T)
end

"""
    properties(A::AnnotatedArray)

Return the object storing the properties associated with the annotated array `A`.

"""
properties(A::AnnotatedArray) = Base.getfield(A, :props)

"""
    properties(T::Type{<:AnnotatedArray})

Return the type of the property object associated with an annotated array of type `T`.

"""
properties(::Type{<:AnnotatedArray{T,N,P}}) where {T,N,P} = P

# Implement syntax `obj.key` for annotated arrays.
function Base.getproperty(A::AbstractAnnotatedArray{<:Any,<:Any,<:SymbolicProperties},
                          key::Symbol)
    return getindex(properties(A), key)
end

# Implement syntax `obj.key = val` for annotated arrays.
function Base.setproperty!(A::AbstractAnnotatedArray{<:Any,<:Any,<:SymbolicProperties},
                           key::Symbol, val)
    setindex!(properties(A), val, key)
    return A
end

function Base.propertynames(A::AbstractAnnotatedArray, private::Bool=false)
    return keys(properties(A))
end

# Method `hasproperty` was introduced in Julia 1.2.
if isdefined(Base, :hasproperty)
    Base.hasproperty(A::AbstractAnnotatedArray, key::Symbol) = hasproperty(properties(A), key)
end

# Substitute for `keytype` and `valtype` that accpet a named tuple or a pair as agrument.
keytype(::Type{T}) where {T<:NamedTuple} = Symbol
valtype(::Type{T}) where {T<:NamedTuple} = eltype(T)
keytype(::Type{Pair{K,V}}) where {K,V} = K
valtype(::Type{Pair{K,V}}) where {K,V} = V
for f in (:keytype, :valtype)
    @eval begin
        $f(x) = Base.$f(x) # fallback to Base method
        $f(x::NamedTuple) = $f(typeof(x))
        $f(x::Pair) = $f(typeof(x))
    end
end

# Substitute for `getkey` that accepts a named tuple as argument.
getkey(A::NamedTuple, key::Symbol, def) = haskey(A, key) ? key : def
getkey(A::NamedTuple, key, def) = def
getkey(A, key, def) = Base.getkey(A, key, def) # fallback to Base method

# Extend methods such as `haskey`, `get`, etc. to access collections by key for annotated
# arrays. We assume that methods like `haskey` and `get!` may be implemented for the values
# of the annotated array even though this is generally not the case for arrays.

function Base.haskey(A::AbstractAnnotatedArray, key)
    if is_key(typeof(properties(A)), typeof(key))
        return haskey(properties(A), key)
    else
        return haskey(values(A), key)
    end
end

function Base.getkey(A::AbstractAnnotatedArray, key, def)
    if is_key(typeof(properties(A)), typeof(key))
        return getkey(properties(A), key, def)
    else
        return getkey(values(A), key, def)
    end
end

# The following are needed to avoid ambiguities.
Base.get(A::AbstractAnnotatedArray, i::Integer, def) = get(values(A), i, def)
Base.get(A::AbstractAnnotatedArray, i::CartesianIndex, def) = get(values(A), i, def)

function Base.get(A::AbstractAnnotatedArray, key, def)
    if is_key(typeof(properties(A)), typeof(key))
        return get(properties(A), key, def)
    else
        return get(values(A), key, def)
    end
end

function Base.get(f::Callable, A::AbstractAnnotatedArray, key)
    val = get(A, key, NotFound())
    val isa NotFound ? f() : val
end


function Base.get!(A::AbstractAnnotatedArray, key, def)
    if is_key(typeof(properties(A)), typeof(key))
        return get!(properties(A), key, def)
    else
        return get!(values(A), key, def)
    end
end

function Base.get!(f::Callable, A::AbstractAnnotatedArray, key)
    if is_key(typeof(properties(A)), typeof(key))
        return get!(f, properties(A), key)
    else
        return get!(f, values(A), key)
    end
end

@propagate_inbounds function Base.getindex(A::AbstractAnnotatedArray, key)
    if is_key(typeof(properties(A)), typeof(key))
        return getindex(properties(A), key)
    else
        return getindex(values(A), key)
    end
end

@propagate_inbounds function Base.setindex!(A::AbstractAnnotatedArray, val, key)
    if is_key(typeof(properties(A)), typeof(key))
        return setindex!(properties(A), val, key)
    else
        return setindex!(values(A), val, key)
    end
end

Base.pop!(A::AbstractAnnotatedArray) = pop!(values(A))
Base.pop!(A::AbstractAnnotatedArray, key) = pop!(properties(A), key)
Base.pop!(A::AbstractAnnotatedArray, key, def) = pop!(properties(A), key, def)

function Base.delete!(A::AbstractAnnotatedArray, key)
    delete!(properties(A), key)
    return A
end

Base.pairs(A::AbstractAnnotatedArray) = pairs(values(A))

"""
    ArrayTools.AnnotatedArrays.is_key(typeof(props), typeof(key))

Return whether `key` is possibly a property name in `props`, the object storing the
properties of an annotated array. The result is a *trait* which only depends on the types of
`props` and of `key`.

To avoid confusing property names and array indices, the default implementation is very
conservative: only strings and symbols may be considered as a property name. However, the
method may be specialized for custom types of `props` or of `keytype(props)` for more
flexibility.

"""
is_key(::Type{<:NamedTuple}, ::Type{T}) where {T} = T === Symbol
is_key(::Type{<:AbstractDict{Symbol}}, ::Type{T}) where {T} = T === Symbol
is_key(::Type{<:AbstractDict{<:AbstractString}}, ::Type{T}) where {T} = T <: AbstractString
is_key(::Type{<:AbstractDict{K}}, ::Type{T}) where {K,T} = false

end # module
