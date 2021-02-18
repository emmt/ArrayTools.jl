#
# AnnotatedArrays.jl -
#
# Objects that combine values stored in an array and properties stored in a
# dictionary or a named tuple.
#

module AnnotatedArrays

export
    AnnotatedArray,
    nkeys,
    properties

using ..ArrayTools.PseudoArrays

# Method `hasproperty` was introduced in Julia 1.2.
if isdefined(Base, :hasproperty)
    import Base: hasproperty
end

"""

`AnnotatedArrays.Properties` is the union of allowed types to store the
properties of `AnnotatedArray` instances.  It can be a `NamedTuple`, a
dictionary with symbolic keys or a dictionary indexed by strings.  The former
case provides the fastest access to properties but they are immutable, the
latter case provides the slowest access to properties.

If other types that symbols or strings are allowed, they must be restricted to
avoid ambiguities with array indices (and ranges, etc.).

For symbolic properties names (stored by a `NamedTuple` or an
`AbstractDict{Symbol}`), the properties can be accessed like ordinary fields.

"""
const Properties = Union{AbstractDict,NamedTuple}
const SymbolicProperties = Union{AbstractDict{Symbol},NamedTuple}

"""

`AbstractAnnotatedArray{T,N,P,S}` is the super-type of annotated arrays which
store values of type `T` as in a `N`-dimensional array with indexing style `S`
and properties in an object of type `P` (a dictionary or a named tuple).

Two useful aliases are defined:

- `DynamicallyAnnotatedArray{T,N,K,V,S}` for annotated arrays whose properties
  can be modified because they are stored in a dictionary with key type `K` and
  value type `V`.

- `StaticallyAnnotatedArray{T,N,S}` for annotated arrays which have immutable
  properties (because they are stored in a named tuple).  To build such
  annotated arrays, provide a named tuple as the last argument of the
  `AnnotatedArray` constructor:

  ```julia
  A = AnnotatedArray(arr, (key=val,))
  B = AnnotatedArray{T}(undef, (dim1, dim2, ...), (key1=val1, key2=val2, ...))
  ```

  Do not forget to have a trailing comma if the named tuple has a single
  element (as for `A` in the above example).

The abstract type `AbstractAnnotatedArray{T,N,P,S}` is defined (but not
exported by `AnnotatedArrays`) so that other types than `AnnotatedArray` can be
derived with similar behavior.  For derived types, say `CustomType`, two
methods should be specialized:

```julia
using ArrayTools.AnnotatedArrays
Base.parent(obj::CustomType) -> the array storing the values of obj
AnnotatedArrays.properties(obj::CustomType) -> the object storing the properties of obj
```

"""
abstract type AbstractAnnotatedArray{T,N,P<:Properties,S} <: PseudoArray{T,N,S} end

const DynamicallyAnnotatedArray{T,N,K,V,S} = AbstractAnnotatedArray{T,N,<:AbstractDict{K,V},S}
@doc @doc(AbstractAnnotatedArray) DynamicallyAnnotatedArray

const StaticallyAnnotatedArray{T,N,S} = AbstractAnnotatedArray{T,N,<:NamedTuple,S}
@doc @doc(AbstractAnnotatedArray) StaticallyAnnotatedArray

struct AnnotatedArray{T,N,P<:Properties,
                      A<:AbstractArray{T,N},
                      S} <: AbstractAnnotatedArray{T,N,P,S}
    data::A
    prop::P
    function AnnotatedArray{T,N,P,A,S}(data::A,
                                       prop::P) where {T,N,P<:NamedTuple,
                                                       A<:AbstractArray{T,N},S}
        @assert IndexStyle(data) === S()
        return new{T,N,P,A,S}(data, prop)
    end
    function AnnotatedArray{T,N,P,A,S}(data::A,
                                       prop::P) where {T,N,K,V,
                                                       P<:AbstractDict{K,V},
                                                       A<:AbstractArray{T,N},S}
        @assert IndexStyle(data) === S()
        K === Any && error("key type must be more specialized than `Any`")
        K <: Integer && error("key type must not be an integer")
        K <: CartesianIndex && error("key type must not be a Cartesian index")
        K <: AbstractRange{<:Integer} && error("key type must not be an integer range")
        K <: Colon && error("key type must not be a colon")
        return new{T,N,P,A,S}(data, prop)
    end
end

# As expected by the PseudoArray interface, extend Base.parent() to return the
# array backing the storage of values.
@inline Base.parent(A::AnnotatedArray) = Base.getfield(A, :data)

# Generic outer constructor with given data and properties.
function AnnotatedArray(data::A, prop::P) where {T,N,P<:Properties,
                                                 A<:AbstractArray{T,N}}
    S = typeof(IndexStyle(data))
    return AnnotatedArray{T,N,P,A,S}(data, prop)
end

# Get rid of matching type-parameter when initial array is specified.  The
# pairs of definitions are needed to disentangle ambiguities.
AnnotatedArray{T,N}(data::AbstractArray{T,N}, args...; kwds...) where {T,N} =
    AnnotatedArray(data, args...; kwds...)
AnnotatedArray{T,N}(data::AbstractArray{T,N}; kwds...) where {T,N} =
    AnnotatedArray(data; kwds...)
AnnotatedArray{T}(data::AbstractArray{T}, args...; kwds...) where {T} =
    AnnotatedArray(data, args...; kwds...)
AnnotatedArray{T}(data::AbstractArray{T}; kwds...) where {T} =
    AnnotatedArray(data; kwds...)

# Constructor for properties provided by keywords.
AnnotatedArray(data::AbstractArray; kwds...) =
    AnnotatedArray(data, Dict{Symbol,Any}(kwds...))

# Constructors with initial array and initial properties specified as key-value
# pairs.  When initial properties are specified as key-value pairs, we want to
# have a dictionary whose key type is specialized if possible (for efficiency)
# but avoid having value type specialized (for flexibility).  If one wants
# specific properties key and value types, it is always possible to explicitly
# specify a dictionary.
AnnotatedArray(data::AbstractArray, args::Pair{K,<:Any}...) where {K} =
    AnnotatedArray(data, Dict{K,Any}(args...))
AnnotatedArray(data::AbstractArray, args::Pair{<:AbstractString,<:Any}...) =
    AnnotatedArray(data, Dict{String,Any}(args...))
AnnotatedArray(data::AbstractArray, args::Pair...) =
    AnnotatedArray(data, _initialproperties(Dict(args...)))

_initialproperties(prop::Dict{K,Any}) where {K} = prop
_initialproperties(prop::Dict{K}) where {K} = convert(Dict{K,Any}, prop)

# Constructors that allocate the initial array.  If only keywords are
# specified, the dimensions of the array may be provided by integer arguments.
AnnotatedArray{T,N}(init, dims::Integer...; kwds...) where {T,N} =
    AnnotatedArray{T,N}(init, dims; kwds...)
AnnotatedArray{T}(init, dims::Integer...; kwds...) where {T} =
    AnnotatedArray{T}(init, dims; kwds...)
AnnotatedArray{T,N}(init, dims::NTuple{N,Integer}, args...; kwds...) where {T,N} =
    AnnotatedArray(Array{T,N}(init, dims), args...; kwds...)
AnnotatedArray{T}(init, dims::Tuple{Vararg{Integer}}, args...; kwds...) where {T} =
    AnnotatedArray(Array{T}(init, dims), args...; kwds...)

# The following ones are to avoid common errors.
AnnotatedArray(::UndefInitializer, dims::Integer...; kwds...) =
    throw_missing_type_parameter()
AnnotatedArray(::UndefInitializer, dims::Tuple{Vararg{Integer}}, args...; kwds...) =
    throw_missing_type_parameter()

"""

```julia
properties(A)
```

yields the properties associated with `A`, an instance of `AnnotatedArray`.

"""
@inline properties(A::AnnotatedArray) = Base.getfield(A, :prop)

"""

```julia
propertyname(T, sym)
```

converts symbol `sym` to a suitable key for an instance of an annotated array
of type `T` (a sub-type of `AbstractAnnotatedArray`), throwing an error if this
conversion is not supported.

"""
propertyname(::Type{<:DynamicallyAnnotatedArray{<:Any,<:Any,Symbol}}, sym::Symbol) = sym
propertyname(::Type{<:StaticallyAnnotatedArray}, sym::Symbol) = sym
propertyname(::Type{<:DynamicallyAnnotatedArray{<:Any,<:Any,String}}, sym::Symbol) =
    String(sym)
@noinline propertyname(::Type{T}, sym::Symbol) where {T<:AbstractAnnotatedArray} =
    error(string("converting symbolic key to ", keytype(T), " is not supported"))

# FIXME: The following specialized method extensions are omitted because
#        `obj.key` is much slower than `obj[key]`
#
#propertyname(::Type{<:DynamicallyAnnotatedArray{<:Any,<:Any,String}}, sym::Symbol) =
#    String(sym)
#propertyname(::Type{<:DynamicallyAnnotatedArray{<:Any,<:Any,K}}, sym::Symbol) where {K} =
#    convert(K, sym)


# Extend methods so that syntax `obj.field` can be used.

@inline Base.getproperty(A::StaticallyAnnotatedArray, sym::Symbol) =
    getproperty(properties(A), sym)

@inline Base.getproperty(A::T, sym::Symbol) where {T<:DynamicallyAnnotatedArray} =
    getindex(properties(A), propertyname(T, sym))

Base.setproperty!(A::StaticallyAnnotatedArray, sym::Symbol, val) =
    throw_immutable_properties()

@inline Base.setproperty!(A::T, sym::Symbol, val) where {T<:DynamicallyAnnotatedArray} =
    setindex!(properties(A), val, propertyname(T, sym))

Base.propertynames(A::StaticallyAnnotatedArray, private::Bool=false) =
    propertynames(properties(A), private)

# FIXME: Should return `Tuple(keys(properties(A)))` to conform to the doc. of
#        `propertynames` but this is slower and for most purposes, an iterable
#        is usually needed.
Base.propertynames(A::DynamicallyAnnotatedArray, private::Bool=false) where {T,N,P} =
    keys(properties(A))


#
# Notes:
#
# * keytype() and valtype() do not take NamedTuple as argument.  So to avoid
#   type-piracy, we define our own keytype() and keytype() methods and just
#   override Base.keytype() and Base.valtype() for AbstractAnnotatedArray.
#
# * The following methods to determine the result of valtype() for a named
#   tuple should cover most of the needs. Maybe we should use the same
#   algorithm as the one used by dictionary constructors to determine the value
#   type.
#
keytype(x) = Base.keytype(x) # use Base definition by default

keytype(::NamedTuple) = Symbol
keytype(::Type{<:NamedTuple}) = Symbol

keytype(::T) where {T<:AbstractAnnotatedArray} = keytype(T)
keytype(::Type{<:AbstractAnnotatedArray{<:Any,<:Any,P}}) where {P} = keytype(P)

Base.keytype(::T) where {T<:AbstractAnnotatedArray} = keytype(T)
Base.keytype(::Type{T}) where {T<:AbstractAnnotatedArray} = keytype(T)

valtype(x) = Base.valtype(x) # use Base definition by default

valtype(::T) where {T<:NamedTuple} = valtype(T)
valtype(::Type{<:NamedTuple{<:Any,T}}) where {T} = _eltype(T)
_eltype(::Type{<:Tuple{Vararg{T}}}) where {T} = T
_eltype(::Type{<:Tuple}) = Any

valtype(::T) where {T<:AbstractAnnotatedArray} = valtype(T)
valtype(::Type{<:AbstractAnnotatedArray{<:Any,<:Any,P}}) where {P} = valtype(P)

Base.valtype(::T) where {T<:AbstractAnnotatedArray} = valtype(T)
Base.valtype(::Type{T}) where {T<:AbstractAnnotatedArray} = valtype(T)

"""

```julia
nkeys(A)
```

yields the number of keys in the dictionary `A` or in the dictionary associated
with `A` if it is an instance of `AttributeArray`.

"""
nkeys(A::AbstractAnnotatedArray) = nkeys(properties(A))
nkeys(A::NamedTuple) = length(A)
nkeys(A::AbstractDict) = length(A)

# Extend methods to access an instance of `AttributeArray` like a dictionary.
Base.haskey(A::DynamicallyAnnotatedArray, key) = haskey(properties(A), key)
Base.haskey(A::StaticallyAnnotatedArray, key::Symbol) = haskey(properties(A), key)
Base.haskey(A::StaticallyAnnotatedArray, key) = false

Base.keys(A::AbstractAnnotatedArray) = keys(properties(A))
Base.values(A::AbstractAnnotatedArray) = values(properties(A))

Base.getkey(A::DynamicallyAnnotatedArray, key, def) = getkey(properties(A), key, def)
Base.getkey(A::StaticallyAnnotatedArray, key::Symbol, def) = (haskey(A, key) ? key : def)
Base.getkey(A::StaticallyAnnotatedArray, key, def) = def

Base.get(A::DynamicallyAnnotatedArray, key, def) = get(properties(A), key, def)
Base.get(A::StaticallyAnnotatedArray, key::Symbol, def) = get(properties(A), key, def)
Base.get(A::StaticallyAnnotatedArray, key, def) = def

Base.get!(A::DynamicallyAnnotatedArray, key, def) = get!(properties(A), key, def)
Base.get!(A::StaticallyAnnotatedArray, key, def) = throw_immutable_properties()

# FIXME: Check whether key type must be narrowed in getindex() and setindex!()
#        to avoid conflict with array indexing.
Base.getindex(A::DynamicallyAnnotatedArray{<:Any,<:Any,K}, key::K) where {K} =
    getindex(properties(A), key)
Base.getindex(A::DynamicallyAnnotatedArray{<:Any,<:Any,String}, key::AbstractString) =
    getindex(properties(A), key)
Base.getindex(A::StaticallyAnnotatedArray, key::Symbol) =
    getindex(properties(A), key)

@inline function Base.setindex!(A::DynamicallyAnnotatedArray{<:Any,<:Any,K},
                                val, key::K) where {K}
    properties(A)[key] = val
    A
end

Base.setindex!(A::StaticallyAnnotatedArray, val, key::Symbol) =
    throw_immutable_properties()

Base.delete!(A::DynamicallyAnnotatedArray, key) = begin
    delete!(properties(A), key)
    return A
end
Base.delete!(A::StaticallyAnnotatedArray, key) = throw_immutable_properties()

Base.pop!(A::DynamicallyAnnotatedArray, key) = pop!(properties(A), key)
Base.pop!(A::DynamicallyAnnotatedArray, key, def) = pop!(properties(A), key, def)
Base.pop!(A::StaticallyAnnotatedArray, key) = throw_immutable_properties()
Base.pop!(A::StaticallyAnnotatedArray, key, def) = throw_immutable_properties()

Base.pairs(A::AbstractAnnotatedArray) = pairs(properties(A))

Base.merge(A::AbstractAnnotatedArray, others...) =
    merge(properties(A), others...)
Base.merge(A::AbstractDict, B::AbstractAnnotatedArray, others...) =
    merge(A, properties(B), others...)
Base.merge(combine::Function, A::AbstractAnnotatedArray, others...) =
    merge(combine, properties(A), others...)
Base.merge(combine::Function, A::AbstractDict, B::AbstractAnnotatedArray, others...) =
    merge(combine, A, properties(B), others...)

Base.merge!(A::DynamicallyAnnotatedArray, others...) = begin
    merge!(properties(A), others...)
    return A
end
Base.merge!(A::StaticallyAnnotatedArray, others...) = throw_immutable_properties()
Base.merge!(A::AbstractDict, B::AbstractAnnotatedArray, others...) =
    merge!(A, properties(B), others...)

Base.merge!(combine::Function, A::DynamicallyAnnotatedArray, others...) = begin
    merge!(combine, properties(A), others...)
    return A
end
Base.merge!(combine::Function, A::StaticallyAnnotatedArray, others...) =
    throw_immutable_properties()
Base.merge!(combine::Function, A::AbstractDict, B::AbstractAnnotatedArray, others...) =
    merge!(combine, A, properties(B), others...)

# Errors.
throw_missing_type_parameter() = error("type parameter T must be specified")
throw_immutable_properties() = error("properties are immutable")

end # module
