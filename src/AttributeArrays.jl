#
# AttributeArrays.jl -
#
# Objects that combine values stored in an array and attribues stored in a dictionary.
#
#-------------------------------------------------------------------------------
#
# This file if part of the TAO software (https://github.com/emmt/ArrayTools)
# licensed under the MIT license.
#
# Copyright (C) 2019, Éric Thiébaut.
#

module AttributeArrays

using ..CopycatArrays

export
    AttributeArray,
    attributes,
    nkeys

struct AttributeArray{T,N,K,V,A<:AbstractArray{T,N},S} <: CopycatArray{T,N,S}
    vals::A
    dict::Dict{K,V}
    function AttributeArray{T,N,K,V,A,S}(vals::A,
                                         dict::Dict{K,V}) where {T,N,K,V,A<:AbstractArray{T,N},S}
        @assert IndexStyle(vals) === S()
        return new{T,N,K,V,A,S}(vals, dict)
    end
end

AttributeArray(vals::AbstractArray, pairs::Pair...) =
    AttributeArray(vals, Dict(pairs...))

function AttributeArray(vals::A, dict::Dict{K,V}) where {T,N,K,V,A<:AbstractArray{T,N}}
    S = typeof(IndexStyle(vals))
    K === Any && error("key type must be more specialized than `Any`")
    K <: Integer && error("key type must not be an integer")
    K <: CartesianIndex && error("key type must not be a Cartesian index")
    return AttributeArray{T,N,K,V,A,S}(vals, dict)
end

AttributeArray{T}(init, dims::NTuple{N,Integer}, pairs::Pair{K,V}...) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), Dict{K,V}(pairs...))
AttributeArray{T}(init, dims::NTuple{N,Integer}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), dict)
AttributeArray{T}(vals::AbstractArray{T,N}, pairs::Pair{K,V}...) where {T,N,K,V} =
    AttributeArray(vals, Dict{K,V}(pairs...))
AttributeArray{T}(vals::AbstractArray{T,N}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(vals, dict)

AttributeArray{T,N}(init, dims::Tuple{Vararg{Integer}}, pairs::Pair{K,V}...) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), Dict{K,V}(pairs...))
AttributeArray{T,N}(init, dims::Tuple{Vararg{Integer}}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), dict)
AttributeArray{T,N}(vals::AbstractArray{T,N}, pairs::Pair{K,V}...) where {T,N,K,V} =
    AttributeArray(vals, Dict{K,V}(pairs...))
AttributeArray{T,N}(vals::AbstractArray{T,N}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(vals, dict)

AttributeArray{T,N,K}(init, dims::Tuple{Vararg{Integer}}, pairs::Pair{K,V}...) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), Dict{K,V}(pairs...))
AttributeArray{T,N,K}(init, dims::Tuple{Vararg{Integer}}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), dict)
AttributeArray{T,N,K}(vals::AbstractArray{T,N}, pairs::Pair{K,V}...) where {T,N,K,V} =
    AttributeArray(vals, Dict{K,V}(pairs...))
AttributeArray{T,N,K}(vals::AbstractArray{T,N}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(vals, dict)

AttributeArray{T,N,K,V}(init, dims::Tuple{Vararg{Integer}}, pairs::Pair...) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), Dict{K,V}(pairs...))
AttributeArray{T,N,K,V}(init, dims::Tuple{Vararg{Integer}}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(Array{T,N}(init, dims), dict)
AttributeArray{T,N,K,V}(vals::AbstractArray{T,N}, pairs::Pair...) where {T,N,K,V} =
    AttributeArray(vals, Dict{K,V}(pairs...))
AttributeArray{T,N,K,V}(vals::AbstractArray{T,N}, dict::Dict{K,V}) where {T,N,K,V} =
    AttributeArray(vals, dict)


# Extend parent() method for the CopycatArray interface to work.
@inline Base.parent(A::AttributeArray) = A.vals

"""

```julia
attributes(A)
```

yields the dictionary associated with `A`, an instance of `AttributeArray`.

"""
@inline attributes(A::AttributeArray) = A.dict

"""

```julia
nkeys(A)
```

yields the number of keys in the dictionary associated with `A`, an instance of
`AttributeArray`.

"""
nkeys(A::AttributeArray) = length(attributes(A))

# Extend methods to access an instance of `AttributeArray` like a dictionary.
Base.haskey(A::AttributeArray, key) = haskey(attributes(A), key)
Base.keys(A::AttributeArray) = keys(attributes(A))
Base.values(A::AttributeArray) = values(attributes(A))
Base.getkey(A::AttributeArray, key, def) = getkey(attributes(A), key, def)
Base.get(A::AttributeArray, key, def) = get(attributes(A), key, def)
Base.get!(A::AttributeArray, key, def) = get!(attributes(A), key, def)
Base.getindex(A::AttributeArray{T,N,K,V}, key::K) where {T,N,K,V} =
    getindex(attributes(A), key)
Base.setindex!(A::AttributeArray{T,N,K,V}, val, key::K) where {T,N,K,V} =
    setindex!(attributes(A), val, key)
Base.delete!(A::AttributeArray, key) = delete!(attributes(A), key)
Base.pop!(A::AttributeArray{T,N,K,V}, key::K) where {T,N,K,V} =
    pop!(attributes(A), key)
Base.pop!(A::AttributeArray{T,N,K,V}, key::K, def) where {T,N,K,V} =
    pop!(attributes(A), key, def)
Base.pairs(A::AttributeArray) = pairs(attributes(A))

Base.merge(A::AttributeArray, others...) =
    merge(attributes(A), others...)
Base.merge(A::AbstractDict, B::AttributeArray, others...) =
    merge(A, attributes(B), others...)
Base.merge(combine::Function, A::AttributeArray, others...) =
    merge(combine, attributes(A), others...)
Base.merge(combine::Function, A::AbstractDict, B::AttributeArray, others...) =
    merge(combine, A, attributes(B), others...)

Base.merge!(A::AttributeArray, others...) =
    merge!(attributes(A), others...)
Base.merge!(A::AbstractDict, B::AttributeArray, others...) =
    merge!(A, attributes(B), others...)
Base.merge!(combine::Function, A::AttributeArray, others...) =
    merge!(combine, attributes(A), others...)
Base.merge!(combine::Function, A::AbstractDict, B::AttributeArray, others...) =
    merge!(combine, A, attributes(B), others...)

Base.keytype(A::AttributeArray{T,N,K,V}) where {T,N,K,V} = K
Base.valtype(A::AttributeArray{T,N,K,V}) where {T,N,K,V} = V

end # module
