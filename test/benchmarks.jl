module ArrayTimings

using Printf
using BenchmarkTools
using BenchmarkTools: prettymemory
using ArrayTools
using LinearAlgebra
using LazyAlgebra

#
# The benchmarks show that `all_indices` is as fast as `eachindex`.
#

flatten(A::AbstractArray{T}) where {T} =
    copyto!(Vector{T}(undef, length(A)), A)

function prettytime(ns)
    a = abs(ns)
    if a < 1e3
        scale, units = 1e0, "ns"
    elseif a < 1e6
        scale, units = 1e-3, "μs"
    elseif a < 1e9
        scale, units = 1e-6, "ms"
    else
        scale, units = 1e-9, "s "
    end
    return string(@sprintf("%8.3f", scale*ns), " ", units)
end

function prettyflops(Gflops)
    a = abs(Gflops)
    if true || a >= 1
        scale, units = 1e0, "Gflops"
    else
        scale, units = 1e3, "Mflops"
    end
    return string(@sprintf("%8.3f", scale*Gflops), " ", units)
end

function prt(str::AbstractString,
             nops::Integer,
             trial::BenchmarkTools.Trial; kwds...)
    prt(stdout, str, nops, trial; kwds...)
end

function prt(io::IO,
             str::AbstractString,
             nops::Integer,
             trial::BenchmarkTools.Trial;
             width::Integer = 0, mark::Bool=false)
    best = minimum(trial)
    print(io, str, repeat(" ", max(1, width - length(str))),
          prettytime(time(best)), "  ",
          prettyflops(nops/time(best)),
          " (", allocs(best), " allocations: ",
          prettymemory(memory(best)), ")")
    mark && print(io,  " [*]")
    print(io,  "\n")
end

#------------------------------------------------------------------------------

sum1(A::AbstractArray) = begin
    s = zero(eltype(A))
    for v in A
        s += v
    end
    return s
end

sum2(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds for v in A
        s += v
    end
    return s
end

sum3(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds @simd for v in A
        s += v
    end
    return s
end

sum4(A::AbstractArray) = begin
    s = zero(eltype(A))
    for i in eachindex(A)
        s += A[i]
    end
    return s
end

sum5(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds for i in eachindex(A)
        s += A[i]
    end
    return s
end

sum6(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds @simd for i in eachindex(A)
        s += A[i]
    end
    return s
end

sum7(A::AbstractArray) = begin
    s = zero(eltype(A))
    for i in all_indices(A)
        s += A[i]
    end
    return s
end

sum8(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds for i in all_indices(A)
        s += A[i]
    end
    return s
end

sum9(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds @simd for i in all_indices(A)
        s += A[i]
    end
    return s
end

#------------------------------------------------------------------------------

function dot4(x::AbstractArray{T,N},
              y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    for i in eachindex(x, y)
        s += x[i]*y[i]
    end
    return s
end

function dot5(x::AbstractArray{T,N},
              y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    @inbounds for i in eachindex(x, y)
        s += x[i]*y[i]
    end
    return s
end

function dot6(x::AbstractArray{T,N},
              y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    @inbounds @simd for i in eachindex(x, y)
        s += x[i]*y[i]
    end
    return s
end

function dot7(x::AbstractArray{T,N},
              y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    for i in all_indices(x, y)
        s += x[i]*y[i]
    end
    return s
end

function dot8(x::AbstractArray{T,N},
              y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    @inbounds for i in all_indices(x, y)
        s += x[i]*y[i]
    end
    return s
end

function dot9(x::AbstractArray{T,N},
              y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    @inbounds @simd for i in all_indices(x, y)
        s += x[i]*y[i]
    end
    return s
end

function unsafe_dot(x::AbstractArray{T,N},
                    y::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    s = zero(T)
    @inbounds @simd for i in eachindex(x)
        s += x[i]*y[i]
    end
    return s
end

#------------------------------------------------------------------------------

norm2_1(A::AbstractArray) = begin
    s = zero(eltype(A))
    for v in A
        s += abs2(v)
    end
    return sqrt(s)
end

norm2_2(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds for v in A
        s += abs2(v)
    end
    return sqrt(s)
end

norm2_3(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds @simd for v in A
        s += abs2(v)
    end
    return sqrt(s)
end

norm2_4(A::AbstractArray) = begin
    s = zero(eltype(A))
    for i in eachindex(A)
        s += abs2(A[i])
    end
    return sqrt(s)
end

norm2_5(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds for i in eachindex(A)
        s += abs2(A[i])
    end
    return sqrt(s)
end

norm2_6(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds @simd for i in eachindex(A)
        s += abs2(A[i])
    end
    return sqrt(s)
end

norm2_7(A::AbstractArray) = begin
    s = zero(eltype(A))
    for i in all_indices(A)
        s +=abs2(A[i])
    end
    return sqrt(s)
end

norm2_8(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds for i in all_indices(A)
        s += abs2(A[i])
    end
    return sqrt(s)
end

norm2_9(A::AbstractArray) = begin
    s = zero(eltype(A))
    @inbounds @simd for i in all_indices(A)
        s += abs2(A[i])
    end
    return sqrt(s)
end

norm2_9(::Type{T}, A::AbstractArray) where {T<:AbstractFloat} = begin
    local s::T = zero(T)
    @inbounds @simd for i in all_indices(A)
        s += abs2(A[i])
    end
    return sqrt(s)
end
#------------------------------------------------------------------------------

function test_sum()
    dims = (100, 100)
    inds = (2:99, 2:99)
    for T in (Float32, Float64)
        A = rand(T, dims)
        V = view(A, inds...)
        X = rand(T, dims)
        Y = rand(T, dims)
        W = view(X, inds...)
        Z = view(Y, inds...)
        a = flatten(A)
        x = flatten(X)
        y = flatten(Y)

        println("\nTest SUM on a regular array (T=$T):")
        n = length(A)  # number of operations
        prt("  sum(a) ", n, @benchmark(sum($a)))
        prt("  sum(A) ", n, @benchmark(sum($A)))
        prt("  sum1(A)", n, @benchmark(sum1($A)))
        prt("  sum2(A)", n, @benchmark(sum2($A)))
        prt("  sum3(A)", n, @benchmark(sum3($A)))
        prt("  sum4(A)", n, @benchmark(sum4($A)))
        prt("  sum5(A)", n, @benchmark(sum5($A)))
        prt("  sum6(A)", n, @benchmark(sum6($A)))
        prt("  sum7(A)", n, @benchmark(sum7($A)))
        prt("  sum8(A)", n, @benchmark(sum8($A)))
        prt("  sum9(a)", n, @benchmark(sum9($a)); mark=true)
        prt("  sum9(A)", n, @benchmark(sum9($A)); mark=true)

        println("\nTest SUM on a sub-array (T=$T):")
        n = length(V)  # number of operations
        prt("  sum ", n, @benchmark(sum($V)))
        prt("  sum1", n, @benchmark(sum1($V)))
        prt("  sum2", n, @benchmark(sum2($V)))
        prt("  sum3", n, @benchmark(sum3($V)))
        prt("  sum4", n, @benchmark(sum4($V)))
        prt("  sum5", n, @benchmark(sum5($V)))
        prt("  sum6", n, @benchmark(sum6($V)))
        prt("  sum7", n, @benchmark(sum7($V)))
        prt("  sum8", n, @benchmark(sum8($V)))
        prt("  sum9", n, @benchmark(sum9($V)); mark=true)
    end
end

function test_dot()
    dims = (100, 100)
    inds = (2:99, 2:99)
    for T in (Float32, Float64)
        A = rand(T, dims)
        V = view(A, inds...)
        X = rand(T, dims)
        Y = rand(T, dims)
        W = view(X, inds...)
        Z = view(Y, inds...)
        a = flatten(A)
        x = flatten(X)
        y = flatten(Y)

        println("\nTest DOT on regular arrays (T=$T):")
        n = 2*length(X)  # number of operations
        prt("  x'*y      ", n, @benchmark($x'*$y))
        prt("  dot(x,y)  ", n, @benchmark(dot($x,$y)))
        prt("  dot4      ", n, @benchmark(dot4($X,$Y)))
        prt("  dot5      ", n, @benchmark(dot5($X,$Y)))
        prt("  dot6      ", n, @benchmark(dot6($X,$Y)))
        prt("  dot7      ", n, @benchmark(dot7($X,$Y)))
        prt("  dot8      ", n, @benchmark(dot8($X,$Y)))
        prt("  dot9      ", n, @benchmark(dot9($X,$Y)); mark=true)
        prt("  unsafe_dot", n, @benchmark(unsafe_dot($X,$Y)))
        prt("  vdot      ", n, @benchmark(vdot($X,$Y)))

        println("\nTest DOT on sub-arrays (T=$T):")
        n = 2*length(W)  # number of operations
        prt("  dot4      ", n, @benchmark(dot4($W,$Z)))
        prt("  dot5      ", n, @benchmark(dot5($W,$Z)))
        prt("  dot6      ", n, @benchmark(dot6($W,$Z)))
        prt("  dot7      ", n, @benchmark(dot7($W,$Z)))
        prt("  dot8      ", n, @benchmark(dot8($W,$Z)))
        prt("  dot9      ", n, @benchmark(dot9($W,$Z)); mark=true)
        prt("  unsafe_dot", n, @benchmark(unsafe_dot($W,$Z)))
        prt("  vdot      ", n, @benchmark(vdot($W,$Z)))
    end
end

function test_norm()
    dims = (100, 100)
    inds = (2:99, 2:99)
    for T in (Float32, Float64)
        A = rand(T, dims)
        V = view(A, inds...)
        X = rand(T, dims)
        Y = rand(T, dims)
        W = view(X, inds...)
        Z = view(Y, inds...)
        a = flatten(A)
        x = flatten(X)
        y = flatten(Y)

        println("\nTest NORM2 on a regular array (T=$T):")
        n = 2*length(A)  # number of operations
        prt("  norm(a,2)         ", n, @benchmark(norm($a,2)))
        prt("  norm(A,2)         ", n, @benchmark(norm($A,2)))
        prt("  norm2_1(A)        ", n, @benchmark(norm2_1($A)))
        prt("  norm2_2(A)        ", n, @benchmark(norm2_2($A)))
        prt("  norm2_3(A)        ", n, @benchmark(norm2_3($A)))
        prt("  norm2_4(A)        ", n, @benchmark(norm2_4($A)))
        prt("  norm2_5(A)        ", n, @benchmark(norm2_5($A)))
        prt("  norm2_6(A)        ", n, @benchmark(norm2_6($A)))
        prt("  norm2_7(A)        ", n, @benchmark(norm2_7($A)))
        prt("  norm2_8(A)        ", n, @benchmark(norm2_8($A)))
        prt("  norm2_9(a)        ", n, @benchmark(norm2_9($a)); mark=true)
        prt("  norm2_9(A)        ", n, @benchmark(norm2_9($A)); mark=true)
        prt("  vnorm2(A)         ", n, @benchmark(vnorm2($A)))
        prt("  norm2_9(Float32,A)", n, @benchmark(norm2_9(Float32,$A)); mark=true)
        prt("  vnorm2(Float32,A) ", n, @benchmark(vnorm2(Float32,$A)))
        prt("  norm2_9(Float64,A)", n, @benchmark(norm2_9(Float64,$A)); mark=true)
        prt("  vnorm2(Float64,A) ", n, @benchmark(vnorm2(Float64,$A)))

        println("\nTest NORM2 on a sub-array (T=$T):")
        n = 2*length(V)  # number of operations
        prt("  norm(V,2)         ", n, @benchmark(norm($V,2)))
        prt("  norm2_1(V)        ", n, @benchmark(norm2_1($V)))
        prt("  norm2_2(V)        ", n, @benchmark(norm2_2($V)))
        prt("  norm2_3(V)        ", n, @benchmark(norm2_3($V)))
        prt("  norm2_4(V)        ", n, @benchmark(norm2_4($V)))
        prt("  norm2_5(V)        ", n, @benchmark(norm2_5($V)))
        prt("  norm2_6(V)        ", n, @benchmark(norm2_6($V)))
        prt("  norm2_7(V)        ", n, @benchmark(norm2_7($V)))
        prt("  norm2_8(V)        ", n, @benchmark(norm2_8($V)))
        prt("  norm2_9(V)        ", n, @benchmark(norm2_9($V)); mark=true)
        prt("  vnorm2(V)         ", n, @benchmark(vnorm2($V)))
        prt("  norm2_9(Float32,V)", n, @benchmark(norm2_9(Float32,$V)); mark=true)
        prt("  vnorm2(Float32,V) ", n, @benchmark(vnorm2(Float32,$V)))
        prt("  norm2_9(Float64,V)", n, @benchmark(norm2_9(Float64,$V)); mark=true)
        prt("  vnorm2(Float64,V) ", n, @benchmark(vnorm2(Float64,$V)))
    end
end

function test_update()
    dims = (100, 100)
    inds = (2:99, 2:99)
    for T in (Float32, Float64)
        X = rand(T, dims)
        Y = similar(X)
        U = view(X, inds...)
        V = view(Y, inds...)
        α = T(0.01)

        println("\nTest UPDATE on regular arrays (T=$T):")
        n = 2*length(X)  # number of operations
        vfill!(Y,0)
        prt("  Y .= Y .+ α.*X   ", n, @benchmark($Y .+= $α.*$X))
        vfill!(Y,0)
        prt("  vupdate!(Y,α,X)  ", n, @benchmark(vupdate!($Y,$α,$X)))

        println("\nTest UPDATE on sub-arrays (T=$T):")
        n = 2*length(U)  # number of operations
        vfill!(V,0)
        prt("  V .= V .+ α.*V   ", n, @benchmark($V .+= $α.*$U))
        vfill!(V,0)
        prt("  vupdate!(V,α,U)  ", n, @benchmark(vupdate!($V,$α,$U)))
    end
end

function test_combine()
    dims = (100, 100)
    inds = (2:99, 2:99)
    for T in (Float32, Float64)
        X = rand(T, dims)
        Y = rand(T, dims)
        Z = similar(X)
        U = view(X, inds...)
        V = view(Y, inds...)
        W = view(Z, inds...)
        α = T(0.3)
        β = T(π)

        println("\nTest COMBINE on regular arrays (T=$T):")
        n = 3*length(X)  # number of operations
        prt("  α*X + β*Y             ", n, @benchmark($α*$X + $β*$Y))
        prt("  vcombine(α,X,β,Y)     ", n, @benchmark(vcombine($α,$X,$β,$Y)))
        prt("  vcombine!(Z,α,X,β,Y)  ", n, @benchmark(vcombine!($Z,$α,$X,$β,$Y)); mark=true)

        println("\nTest COMBINE on sub-arrays (T=$T):")
        n = 3*length(U)  # number of operations
        prt("  α*U + β*V             ", n, @benchmark($α*$U + $β*$V))
        prt("  vcombine(α,U,β,V)     ", n, @benchmark(vcombine($α,$U,$β,$V)))
        prt("  vcombine!(W,α,U,β,V)  ", n, @benchmark(vcombine!($W,$α,$U,$β,$V)); mark=true)
    end
end

end # module

if !isinteractive()
    println("\nIn these tests, the asterisk [*] marks the version ",
            "which should be among the fastest ones and which should ",
            "require no additional memory.")
    ArrayTimings.test_sum()
    ArrayTimings.test_dot()
    ArrayTimings.test_norm()
    ArrayTimings.test_update()
    ArrayTimings.test_combine()
end

nothing
