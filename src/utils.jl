#
# utils.jl --
#
# General purpose methods.
#

"""
    to_type(T, x)

yields `x` converted to type `T`, the result is asserted to be of type `T`.

"""
to_type(::Type{T}, x::T) where {T} = x
to_type(::Type{T}, x::Any) where {T} = convert(T, x)::T

"""
    to_int(x)

converts `x` to a similar object whose values are all of type `Int`.  Argument
`x` can be a scalar, an array, a tuple or a range.  Argument `x` is returned if
already of the correct type.

"""
to_int(x::Int) = x
to_int(x::Integer) = Int(x)
to_int(I::AbstractUnitRange{Int}) = I
to_int(I::AbstractUnitRange{<:Integer}) = to_int(first(I)):to_int(last(I))
to_int(I::StepRange{Int,Int}) = I
to_int(I::StepRange{<:Integer,<:Integer}) =
    to_int(first(I)):to_int(step(I)):to_int(last(I))
to_int(A::AbstractArray{Int,N}) where {N} = A
to_int(A::AbstractArray{<:Integer,N}) where {N} = map(to_int, A)
to_int(T::NTuple{N,Int}) where {N} = T
to_int(T::NTuple{N,Integer}) where {N} = map(to_int, T)

"""
    promote_eltype(A, B, ...) -> T

yields an element type `T` resulting from calling `promote_type` onto the
element types of all the arguments which must be arrays or array types.

"""
promote_eltype(args::Union{AbstractArray,Type{<:AbstractArray}}...) =
    promote_type(map(eltype, args)...)

"""
    all_match(val, f, args...) -> bool

yields as soon as possible (short-circuit) whether `f(arg) == val` for each
argument `arg` in `args...`.  The returned value is `true` if there are no
arguments after `f`.

"""
all_match(val, f::Function) = true
all_match(val, f::Function, A) = f(A) == val
@inline all_match(val, f::Function, A, B...) =
    all_match(val, f, A) && all_match(val, f::Function, B...)

"""
    allof(f, args...) -> Bool

checks whether predicate function `f` returns `true` for all arguments in
`args...`, returning `false` as soon as possible (short-circuiting).

    allof(args...) -> Bool

checks whether all arguments `args...` are `true`, returning `false` as soon as
possible (short-circuiting).  Arguments can be booleans or arrays of booleans.
The latter are considered as `true` if all their elements are `true` and are
considered as `false` otherwise (if any of their elements are `false`).
Arguments can also be iterables to check whether all their values are `true`.
An empty iterable is considered as `true`.

This method can be much faster than `all(f, args)` or `all(args)` because its
result may be determined at compile time.  However, `missing` values are not
considered as special.

See also `all`, [`anyof`](@ref), [`noneof`](@ref).

"""
allof(f::Function, a) = f(a)::Bool
@inline allof(f::Function, a, b...) = allof(f, a) && allof(f, b...)
@inline allof(a, b...) = allof(a) && allof(b...)
allof(a::Bool) = a
function allof(a::AbstractArray{Bool})
    @inbounds for i in eachindex(a)
        a[i] || return false
    end
    return true
end
function allof(itr)
    for val in itr
        allof(val) || return false
    end
    return true
end

"""
    anyof(f, args...) -> bool

checks whether predicate function `f` returns `true` for any argument
`args...`, returning `true` as soon as possible (short-circuiting).

    anyof(args...) -> bool

checks whether all arguments `args...` are `true`, returning `false` as soon as
possible (short-circuiting).  Arguments can be booleans or arrays of booleans.
The latter are considered as `true` if any of their elements are `true` and are
considered as `false` otherwise (if all their elements are `false`).  Arguments
can also be iterables to check whether any of their values are `true`.  An
empty iterable is considered as `false`.

This method can be much faster than `any(f, args)` or `any(args)` because its
result may be determined at compile time.  However, `missing` values are not
considered as special.

See also `any`, [`allof`](@ref), [`noneof`](@ref).

"""
anyof(f::Function, a) = f(a)::Bool
@inline anyof(f::Function, a, b...) = anyof(f, a) || anyof(f, b...)
@inline anyof(a, b...) = anyof(a) || anyof(b...)
anyof(a::Bool) = a
function anyof(a::AbstractArray{Bool})
    @inbounds for i in eachindex(a)
        a[i] && return true
    end
    return false
end
function anyof(itr)
    for val in itr
        anyof(val) && return true
    end
    return false
end

"""
    noneof(f, args...) -> bool

checks whether predicate `f` returns `false` for all argument `args...`, while

    noneof(args...) -> bool

checks whether all argument `args...` are false.

See also `any`, [`allof`](@ref), [`noneof`](@ref).

"""
@inline noneof(args...) = ! anyof(args...)

"""
    reversemap(f, args)

applies the function `f` to arguments `args` in reverse order and return the
result.  For now, the arguments `args` must be in the form of a simple tuple
and the result is the tuple: `(f(args[end]),f(args[end-1]),...,f(args[1])`.

Also see: `map`, `ntuple`.

"""
reversemap(f::Function, args::NTuple{N,Any}) where {N} =
    ntuple(i -> f(args[(N + 1) - i]), Val{N}())

"""
    strictmap!(dst, f, src) -> dst

does `dst[i] = f(src[i])` for all indices `i` and returns `dst`.  Arguments
`dst` and `src` must have the same axes.

Except for the strict condition on the axes, this method is similar to
`map!(f,dst,src)`.

"""
function strictmap!(dst::AbstractArray{<:Any,N}, f,
                    src::AbstractArray{<:Any,N}) where {N}
    @inbounds @simd for i in all_indices(dst, src)
        dst[i] = f(src[i])
    end
    return dst
end
