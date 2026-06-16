module ArrayToolsTests

using Test, Random
using ArrayTools, ArrayTools.AnnotatedArrays, ArrayTools.PseudoArrays
using Base: @propagate_inbounds
import Base: getindex, setindex!, checkbounds

function samevalues(A::AbstractArray, B::AbstractArray)
    @assert has_standard_indexing(A, B)
    @assert size(A) == size(B)
    for i ∈ eachindex(A, B)
        A[i] != B[i] && return false
    end
    return true
end

function maxabsdif(A::AbstractArray, B::AbstractArray)
    @assert has_standard_indexing(A, B)
    @assert size(A) == size(B)
    T = promote_type(eltype(A), eltype(B))
    res = zero(T)
    for i ∈ eachindex(A, B)
        res = max(res, abs(A[i] - B[i]))
    end
    return res
end

slice(A::AbstractArray{T,N}, i::Integer) where {T,N} =
    A[colons(Val(N-1))..., i]

generate(::Type{T}, dims::Integer...) where {T} = generate(T, dims)
generate(::Type{T}, dims::NTuple{N,Integer}) where {T,N} =
    generate!(Array{T,N}(undef, dims))
generate!(A::AbstractArray{T,N}) where {T,N} = begin
    k = 0
    @inbounds for i in eachindex(A)
        k += 1
        A[i] = k
    end
    return A
end

const T = Float32
const dims = (3, 4, 5, 6)
const N = length(dims)
A = generate(T, dims)
Vf = view(A, :, :, :, 3:4)   # a flat sub-array
Va = view(A, 2:3, :, :, 3:4) # a non-flat sub-array
V = view(A, :, 2:3, :, :)
S = 1:2:70  # StepRange
U = 3:50    # UnitRange
atol = 1e-6

@testset "Miscellaneous " begin
    # Promotion of element types.
    for (T1,T2,T3) in ((Float32,Float64,Int), (Int16,Int32,Int64))
        @test promote_eltype() === promote_type()
        @test promote_eltype(zeros(T1,1)) == promote_type(T1)
        @test promote_eltype(AbstractVector{T2}) == promote_type(T2)
        @test promote_eltype(zeros(T2,4), zeros(T3,2,3)) == promote_type(T2,T3)
        @test promote_eltype(zeros(T2,2), Array{T3}) == promote_type(T2,T3)
        @test promote_eltype(zeros(T1,3), zeros(T2,5), AbstractVector{T3}) ==
            promote_type(T1,T2,T3)
        @test promote_eltype(zeros(T1,3), zeros(T2,5), ones(T3,2)) ==
            promote_type(T1,T2,T3)
        @test promote_eltype(DenseMatrix{T1}, zeros(T2,5), AbstractVector{T3}) ==
            promote_type(T1,T2,T3)
    end
    # Conversions.
    @test_deprecated to_type(Int32, 3)
    @test as(Int32, 3) === Int32(3)
    @test as(Float32, 3) === Float32(3)
    @test as(Vector{Float32}, [1,2,3]) == Float32[1,2,3]
    @test as(Vector{Float32}, [1,2,3]) isa Vector{Float32}
    @test to_int(2) === 2
    @test to_int(0x03) === 3
    @test to_int(Base.OneTo(5)) === Base.OneTo(5)
    @test to_int(Base.OneTo{Int16}(5)) === Base.OneTo(5)
    @test to_int(-1:5) === -1:5
    @test to_int(0x00:0x5) === 0:5
    @test to_int(-1:2:7) === -1:2:7
    @test to_int(0x00:0x02:0x09) === 0:2:8
    @test to_int(Int8(9):Int8(-3):Int8(-1)) === 9:-3:0
    let A = [1, -2, 11]
        @test to_int(A) === A
        @test to_int(Int16.(A)) == A
    end
    let A = (1, -2, 11)
        @test to_int(A) === A
        @test to_int(Int16.(A)) == A
    end

    # Dimensions.
    @test_deprecated to_size(3) == (3,)
    @test standard_size(A) === size(A)
    @test same_standard_size(A) === size(A)
    @test same_standard_size(A, similar(A)) === size(A)
    @test_throws DimensionMismatch same_standard_size(A, V)
    @test same_size(A) === size(A)
    @test same_size(A, similar(A)) === size(A)
    @test_throws DimensionMismatch same_size(A, V)
    @test check_size(()) == 1
    @test check_size((1,0,2)) == 0
    @test check_size((1,2,3)) == 6
    @test check_size(0xff) == 255
    @test_throws ArgumentError check_size((1,0,-1))
    @test_throws ArgumentError check_size(-1)
    # Axis and axes.
    @test_deprecated to_axis(1:3) === 1:3
    @test_deprecated to_axes(1:3) === (1:3,)
    @test_deprecated to_axes(2:3,7) === (2:3,Base.OneTo(7))
    #
    # Tests for `allof`, `anyof` and `noneof`.
    #
    W = [true, true, true]
    X = [true, false]
    Y = [false, false]
    for i in randperm(length(X)) # prevent compilation-time optimization
        @test anyof(X[i]) == X[i]
    end
    @test allof(true) == true
    @test anyof(true) == true
    @test noneof(true) == false
    @test allof(false) == false
    @test anyof(false) == false
    @test noneof(false) == true
    @test allof(W) == true
    @test anyof(W) == true
    @test noneof(W) == false
    @test allof(Tuple(W)) == allof(W)
    @test anyof(Tuple(W)) == anyof(W)
    @test noneof(Tuple(W)) == noneof(W)
    @test allof(X) == false
    @test anyof(X) == true
    @test noneof(X) == false
    @test allof(Tuple(X)) == allof(X)
    @test anyof(Tuple(X)) == anyof(X)
    @test noneof(Tuple(X)) == noneof(X)
    @test allof(Y) == false
    @test anyof(Y) == false
    @test noneof(Y) == true
    @test allof(Tuple(Y)) == allof(Y)
    @test anyof(Tuple(Y)) == anyof(Y)
    @test noneof(Tuple(Y)) == noneof(Y)
    for (x, y) in ((W,W),(X,X),(Y,Y),(W,X),(X,Y),(Y,W))
        @test anyof(x, y) == (anyof(x) || anyof(y))
        @test allof(x, y) == (allof(x) && allof(y))
        @test noneof(x, y) == (noneof(x) && noneof(y))
    end
    for x in (W, X, Y)
        @test allof(x) == (allof(minimum, x, x) && allof(maximum, x))
        @test noneof(x) == (noneof(minimum, x, x) && noneof(maximum, x))
        @test anyof(x) == (anyof(minimum, x, x) || anyof(maximum, x))
    end
    #
    # Tests for `cartesian_indices`.
    #
    @test cartesian_indices(A) === CartesianIndices(axes(A))
    @test cartesian_indices(A) === cartesian_indices(size(A))
    @test cartesian_indices(A) === cartesian_indices(axes(A))
    @test cartesian_indices(A) === cartesian_indices(cartesian_indices(A))
    I1 = CartesianIndex(1,2,3)
    I2 = CartesianIndex(5,7,9)
    @test cartesian_indices(I1,I2) ===
        CartesianIndices(([I1[k]:I2[k] for k in 1:length(I1)]...,))
    @test cartesian_indices(3:4) === CartesianIndices((3:4,))
    @test cartesian_indices(3:4,2:5) === CartesianIndices((3:4,2:5))
    @test cartesian_indices((3:4,2:5)) === CartesianIndices((3:4,2:5))
    #
    # Tests for `axis_limits`.
    #
    @test axis_limits(Base.OneTo(7)) === (1,7)
    @test axis_limits(7:16) === (7,16)
    @test axis_limits(Int16(7):Int16(1):Int16(16)) === (7,16)
    @test axis_limits(16:-1:7) === (7,16)
    @test axis_limits(7:-1:16) === (8,7)
    @test_throws ArgumentError axis_limits(7:3:16)
    #
    # Tests for `all_match`.
    #
    @test all_match(identity, nothing)
    @test all_match(length, 1, π)
    @test all_match(length, 1, π, 2, "e")
    @test all_match(NaN, :w, 1, π) do _; NaN; end
    @test_deprecated all_match(1, length, π)
    #
    # Tests for `same_axes`.
    #
    let A = Array{Char}(undef, (2,3)),
        B = Array{Int}(undef, 2, 3),
        C = Array{Float32}(undef, (2,3,4))
        @test same_axes(A) === axes(A)
        @test same_axes(B) === axes(B)
        @test same_axes(A,B) == axes(A)
        @test same_axes(A,B) == axes(B)
        @test_throws DimensionMismatch same_axes(A,B,C)
    end
    #
    # Tests for `all_indices`.
    #
    B = rand(T, dims[1:2])
    C = rand(T, dims[1:end-2]..., dims[end], dims[end-1])
    X = rand(T, 6)
    @test all_indices(X) === eachindex(X)
    @test all_indices(A) === eachindex(A)
    @test all_indices(Va) === eachindex(Va)
    @test_throws DimensionMismatch all_indices(A,Va)
    @test eachindex(A,C) === eachindex(A)
    @test_throws DimensionMismatch all_indices(A,C)
    @test all_indices(A, rand(T, dims), rand(T, dims)) === eachindex(A)
    Y = rand(T, dims[1], dims[2]+2, dims[3:end]...)
    Z = view(Y, :, 2:dims[2]+1, colons(length(dims)-2)...)
    @test IndexStyle(Z) === IndexCartesian()
    @test all_indices(A, rand(T, dims), Z) === eachindex(IndexCartesian(), A)
    #
    # Other stuff.
    #
    @test reversemap(x -> x^2, dims) === reverse(map(x -> x^2, dims))
    @test_throws ArgumentError ArrayTools.throw_non_standard_indexing()
    let A = rand(Float32, (2,3,4)),
        B = Array{Float32}(undef, size(A)),
        C = Array{Float32}(undef, map(d -> d + 1, size(A)))
        @test strictmap!(-, B, A) === B
        @test B == map(-, A)
        @test strictmap!(B, -, A) === B
        @test B == map(-, A)
        @test_throws DimensionMismatch strictmap!(-, C, A)
        @test_throws ArgumentError strictmap!(A, B, C)
    end
end

#
# Tests for rubber indices.
#
@testset "Rubber Indices" begin
    I1 = CartesianIndex(2)
    I2 = CartesianIndex(2,3)
    I3 = CartesianIndex(2,3,4)
    I4 = CartesianIndex(1,2,3,2)
    I5 = CartesianIndex(1,2,3,2,4)
    @test_throws ArgumentError colons(-1)
    @testset "$d colon(s)" for d ∈ 0:15
        tup = ntuple(x -> Colon(), d)
        @test colons(d) === tup
        @test colons(Val(d)) === tup
    end
    for k ∈ 1:dims[end]
        @test samevalues(slice(A, k), A[:,:,:,k])
    end
    @test A[..]           === A
    @test A[..]           == A[:,:,:,:]
    @test A[..,3]         == A[:,:,:,3]
    @test A[2,..]         == A[2,:,:,:]
    @test A[..,2:4,5]     == A[:,:,2:4,5]
    @test A[2:3,..,1,2:4] == A[2:3,:,1,2:4]
    @test A[:,2,..,2:4]   == A[:,2,:,2:4]
    @test A[:,2:3,..]     == A[:,2:3,:,:]
    @test A[..,2:3,:]     == A[:,:,2:3,:]
    @test A[I1,..]        == A[I1,:,:,:]
    @test A[I2,..]        == A[I2,:,:]
    @test A[..,I1]        == A[:,:,:,I1]
    @test A[..,I2]        == A[:,:,I2]
    @test A[I1,..,I2]     == A[I1,:,I2]
    @test A[I2,..,I1]     == A[I2,:,I1]
    @test A[I1,..,I3]     == A[I1,I3]
    @test A[1,2:end,..]   == A[1,2:end,:,:]
    @test A[1,2:end-1,..] == A[1,2:end-1,:,:]
    # The `end` keyword can only appear **before** the rubber index.
    @test_broken A[2,..,3:end] == A[2,:,:,3:end]
    A[1,..] .= 3
    @test all(isequal(3), A[1,:,:,:])
    A[1,..] = A[2,..]
    @test A[1,:,:,:] == A[2,:,:,:]
    A[1,..,2:3] .= A[2,..,3:2:5]
    @test A[1,:,:,2:3] == A[2,:,:,3:2:5]
    # Tests with a 9 dimensional array.
    siz = (3, 4, 2, 3, 4, 2, 3, 4, 2)
    n = length(siz)
    B = generate(Int, siz)
    C = copy(B)
    inds = (2:3, 4, Colon(), CartesianIndex(2),
            CartesianIndex(3), 1:2, 3, 2:3, CartesianIndex(1))
    for k in 0:n, l in 0:n-k
        J1 = inds[1:k]         # leading indices
        J2 = inds[n-l+1:n]     # trailing indices
        J = (J1..., .., J2...,)# with a rubber index
        K = (J1..., colons(n - length(J1) - length(J2))..., J2...,)
        R = B[K...]            # extract region before change
        @test R == B[J...]     # extract values with `getindex`
        B[J...] .*= -1         # change values with `dotview`
        @test R == -B[J...]    # extract values with `getindex`
        B[J...] = R            # restore values with `setindex!`
        @test B == C
    end
end

@testset "Storage       " begin
    B = to_flat_array(Float32, A)
    C = to_flat_array(Float32, Va)
    @test StorageType() === AnyStorage()
    @test StorageType("a") === AnyStorage()
    @test StorageType(A) === FlatStorage()
    @test StorageType(Vf) === FlatStorage()
    @test StorageType(Va) === AnyStorage()
    @test is_flat_array() == false
    @test is_flat_array("a") == false
    @test is_flat_array(S) == false
    @test is_flat_array(A) == true
    @test is_flat_array(Vf) == true
    @test is_flat_array(Va) == false
    @test is_flat_array(A,B,C) == true
    @test is_flat_array(A,B,C,Vf) == true
    @test is_flat_array(A,B,C,Va) == false
    @test to_flat_array(A) === A
    @test to_flat_array(eltype(A), A) === A
    @test pointer(A) != pointer(to_flat_array(Va))
    @test samevalues(S, to_flat_array(S))
    @test eltype(S) === eltype(to_flat_array(S))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(Va, C) ≤ atol
    Da = to_flat_array(Va)
    @test is_flat_array(Da) == true
    @test Da == Va
    @test Da !== Va
    Df = to_flat_array(Vf)
    @test is_flat_array(Df) == true
    @test Df == Vf
    @test Df === Vf
    for n in 1:5
        K = rand(Float64, ntuple(x -> 3, n))
        L = view(K, 2:3, colons(n-1)...)
        M = view(K, colons(n-1)..., 2:3)
        @test StorageType(K) == FlatStorage()
        @test StorageType(L) == (n == 1 ? FlatStorage() : AnyStorage())
        @test StorageType(M) == FlatStorage()
    end
    L = ((nothing,               false),
         ("a",                   false),
         ((),                    false),
         (1,                     false),
         (A,                     true ),
         (view(A, :, :, :, 2:3), true ),
         (view(A, :, :, 2:2, :), false),
         (view(A, :, :, 2, :),   false),
         (view(A, :, :, 2:2, 3), true ),
         (view(A, :, :, 2, 3),   true ),
         (view(A, :, :, 2, 3:3), false))
    for i in randperm(length(L)) # prevent compilation-time optimization
        x, b = L[i]
        @test is_flat_array(x) == b
    end
    @test is_flat_array(A, view(A, :, :, 2, 3), view(A, :, :, 2:2, 3)) == true
    @test is_flat_array(A, view(A, :, :, 2:2, :), view(A, :, :, 2:2, 3)) == false
end

@testset "Indexing      " begin
    # Prior to version 1.1, Base.has_offset_axes does not work for tuples.
    for Q in (A,Va,S,101,Colon(),(),(1,),("e",2,))
        if isa(Q, Tuple) && VERSION < v"1.1"
            @test has_standard_indexing(Q) == true
        else
            @test has_standard_indexing(Q) == !Base.has_offset_axes(Q)
        end
    end
    B = to_fast_array(Float32, A)
    C = to_fast_array(Float32, Va)
    @test has_standard_indexing(A,Va) == (has_standard_indexing(A) &&
                                          has_standard_indexing(Va))
    @test IndexingType(A) === FastIndexing()
    @test IndexingType(Va) === AnyIndexing()
    @test IndexingType("a") === AnyIndexing()
    @test is_fast_array() == false
    @test is_fast_array(S) == true
    @test samevalues(Va, to_fast_array(Va))
    @test samevalues(A, to_fast_array(A))
    @test is_fast_array(A) == true
    @test is_fast_array(A,B,C) == true
    @test is_fast_array(Va) == false
    @test is_fast_array(A,B,C,Va) == false
    @test is_fast_array(to_fast_array(Va)) == true
    @test pointer(A) == pointer(to_fast_array(A))
    @test pointer(A) != pointer(to_fast_array(Va))
    @test samevalues(Va, to_fast_array(Va))
    @test samevalues(A, to_fast_array(A))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(Va, C) ≤ atol

    let w = 0:2, x = 1:3, y = 2:4, z = -1:5
        @test (@assert_same_axes) === nothing
        @test (@assert_same_axes w) === nothing
        @test (@assert_same_axes w x) === nothing
        @test (@assert_same_axes w x y) === nothing
        @test_throws DimensionMismatch (@assert_same_axes w x y z)
    end
end

@testset "Broadcasting  " begin
    # Check broadcasting of dimensions.
    @test_throws DimensionMismatch bcastsize(3, 4)
    for (a,b) in (((), ()),
                  ((Int16(5), Int32(6), Int64(7)), (5,6,7)),)
        @test bcastsize(a) === b
        @test bcastsize((), a) === b
        @test bcastsize(a, ()) === b
        @test bcastsize(a, (a..., 3,)) === (b..., 3,)
        @test bcastsize(a, (a..., 3,), (a..., 1, 2)) === (b..., 3, 2)
        @test bcastsize((a..., 3,), a) === (b..., 3,)
        @test bcastsize((a..., 1,), (a..., 3,)) === (b..., 3,)
        @test bcastsize((a..., 3,), (a..., 1,)) === (b..., 3,)
    end
    # Check `bcastcopy`.
    T = Int32
    dims1 = (1,2,3)
    dims2 = (3,2,1,4)
    dims3 = bcastsize(dims1, dims2)
    A1 = generate(T, dims1)
    A2 = generate(T, dims2)
    @test A1 .+ zeros(T, dims2) == bcastcopy(A1, dims3...)
    @test A2 .+ zeros(T, dims1) == bcastcopy(A2, map(Int16, dims3))
    C1 = bcastcopy(A, eltype(A), size(A))
    @test C1 !== A && C1 == A
    C2 = bcastcopy(A, eltype(A), map(Int16, size(A))...)
    @test C2 !== A && C2 == A
    C4 = bcastcopy(A, size(A))
    @test C4 !== A && C4 == A
    C5 = bcastcopy(A, map(Int16, size(A))...)
    @test C5 !== A && C5 == A
    # Check that `bcastlazy attempt` to yield the same object.
    @test bcastlazy(A, eltype(A), size(A)) === A
    @test bcastlazy(A, eltype(A), map(Int16, size(A))...) === A
    @test bcastlazy(A, size(A)) === A
    @test bcastlazy(A, map(Int16, size(A))...) === A
    B1 = bcastlazy(A, Int32, size(A))
    @test B1 !== A && B1 == A
    B2 = bcastlazy(Float32(1), dims)
    @test eltype(B2) === Float32
    @test size(B2) === dims
    @test B2 == ones(Float32, dims)
    B3 = bcastlazy(1, Float16, dims...)
    @test eltype(B3) === Float16
    @test size(B3) === dims
    @test B3 == ones(Float16, dims)
    # Check slices of arrays given by `bcastcopy`/`bcastlazy`
    X = generate(Float64, (dims[1], 1, dims[3:end]...))
    X1 = bcastlazy(X, dims)
    X2 = bcastcopy(X, dims)
    for i in 1:dims[2]
        @test X[:,1,:,:] == X1[:,i,:,:]
        @test X[:,1,:,:] == X2[:,i,:,:]
    end
end

# UnfinishedArray does not extend Base.parent()
struct UnfinishedArray{T,N,A<:AbstractArray{T,N},S} <: PseudoArray{T,N,S}
    arr::A
    cnt::Int
end
UnfinishedArray(arr::A, cnt::Integer=0) where {T,N,A<:AbstractArray{T,N}} = begin
    S = typeof(IndexStyle(arr))
    return UnfinishedArray{T,N,A,S}(arr,cnt)
end

# DummyArray does extend Base.parent()
struct DummyArray{T,N,A<:AbstractArray{T,N},S} <: PseudoArray{T,N,S}
    arr::A
    cnt::Int
end
DummyArray(arr::A, cnt::Integer=0) where {T,N,A<:AbstractArray{T,N}} = begin
    S = typeof(IndexStyle(arr))
    return DummyArray{T,N,A,S}(arr,cnt)
end
Base.parent(A::DummyArray) = A.arr

@testset "Annotated arrays" begin
    inds = map(n -> Base.OneTo(n), dims)
    T = Float32
    N = length(dims)

    # Forbidden key types.
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{Any,Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{Int32,Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{CartesianIndex,Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{UnitRange{Int},Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{Colon,Any}())

    # Common errors.
    @test_throws Exception AnnotatedArray(undef, dims)
    @test_throws Exception AnnotatedArray(undef, dims...)

    # Annotated arrays with immutable properties.
    props = (a = 1, b = 2, c = 3)
    @test @inferred(AnnotatedArrays.keytype(props)) === Symbol
    @test @inferred(AnnotatedArrays.keytype(typeof(props))) === Symbol
    @test @inferred(AnnotatedArrays.valtype(props)) === Int
    @test @inferred(AnnotatedArrays.valtype(typeof(props))) === Int
    X = @inferred(AnnotatedArray(A, props))
    @test isa(X, AnnotatedArray{T,N})
    # Accessors.
    @test @inferred(values(X)) === A
    @test @inferred(properties(X)) === props
    # Other constructors.
    @test @inferred(AnnotatedArray{T}(A, props)) === X
    @test @inferred(AnnotatedArray{T,N}(A, props)) === X
    @test_throws Exception AnnotatedArray(undef, size(A), props) # NOTE missing element type `T`
    Y = @inferred(AnnotatedArray{T}(undef, size(A), props))
    @test isa(Y, AnnotatedArray{T,N})
    @test @inferred(values(Y)) isa Array{T,N}
    @test @inferred(properties(Y)) === props
    Z = @inferred(AnnotatedArray{T,N}(undef, size(A), props))
    @test isa(Z, AnnotatedArray{T,N})
    @test @inferred(values(Z)) isa Array{T,N}
    @test @inferred(properties(Z)) === props
    # Test `keytype` and `valtype`.
    @test @inferred(keytype(X)) === @inferred(keytype(typeof(X))) === keytype(A)
    @test @inferred(valtype(X)) === @inferred(valtype(typeof(X))) === valtype(A)
    @test @inferred(AnnotatedArrays.keytype(properties(X))) === Symbol
    @test @inferred(AnnotatedArrays.keytype(properties(typeof(X)))) === Symbol
    @test @inferred(AnnotatedArrays.valtype(properties(X))) === Int
    @test @inferred(AnnotatedArrays.valtype(properties(typeof(X)))) === Int
    # Abstract array API.
    @test @inferred(length(X)) == length(values(X)) == length(A)
    @test @inferred(size(X)) === size(A)
    @test @inferred(axes(X)) === axes(A)
    @test @inferred(ndims(X)) === @inferred(ndims(typeof(X))) === ndims(A)
    @test @inferred(eltype(X)) === @inferred(eltype(typeof(X))) === eltype(A)
    @test @inferred(firstindex(X)) === firstindex(A)
    @test @inferred(lastindex(X)) === lastindex(A)
    @test @inferred(CartesianIndices(X)) == CartesianIndices(A)
    @test @inferred(IndexStyle(X)) === IndexStyle(A)
    @test mapreduce(isequal, &, X, A)
    @test @inferred(pairs(X)) === pairs(A)
    @test X == A # only compare values (FIXME?)
    @test @inferred(copyto!(Y, 1:length(Y))) === Y
    @test @inferred(fill!(Z, zero(T))) === Z
    @test Z != Y
    for i in eachindex(Y, Z)
        Z[i] = Y[i]
    end
    @test Z == Y
    # Properties.
    @test @inferred(properties(X)) isa NamedTuple
    @test length(properties(X)) == 3
    let names = @inferred(propertynames(X))
        @test length(names) == length(properties(X)) == 3
        @test (:a ∈ names) && (:b ∈ names) && (:c ∈ names)
    end
    @test X.a === props.a
    @test X.b === props.b
    @test X.c === props.c
    @test_throws Exception X.unknown
    # Indexation by array index or by property key.
    i_start, i_stop = firstindex(A), lastindex(A)
    rng = i_start : min(i_start + 3, i_stop)
    I_start, I_stop = first(CartesianIndices(A)), last(CartesianIndices(A))
    I_step = CartesianIndex(ntuple(i -> i == 1 ? 1 : 0, Val(N)))
    @test X[i_start] == A[i_start]
    @test X[i_stop] == A[i_stop]
    @test X[I_start] == A[I_start]
    @test X[I_stop] == A[I_stop]
    @test X[rng] == A[rng]
    @test (X[:a], X[:b], X[:c]) === (props.a, props.b, props.c)
    @test  haskey(X, :a)
    @test !haskey(X, :x)
    @test_throws MethodError haskey(X, i_start) # NOTE `haskey` is not implemented for arrays
    @test getkey(X, :b, undef) === :b
    @test getkey(X, :x, undef) === undef
    @test_throws MethodError getkey(X, i_start, undef) # NOTE `getkey` is not implemented for arrays
    @test get(X, :c, undef) == X.c
    @test get(X, :x, undef) == undef
    @test get(X, firstindex(X), undef) == first(A)
    @test get(X, firstindex(X) - 1, undef) === undef
    @test get(X, lastindex(X), undef) === last(A)
    @test get(X, lastindex(X) + 1, undef) === undef
    @test get(X, I_start, undef) === first(A)
    @test get(X, I_start - I_step, undef) === undef
    @test get(X, I_stop, undef) === last(A)
    @test get(X, I_stop + I_step, undef) === undef
    @test_throws Exception delete!(X, :a) # NOTE X properties are immutable
    @test_throws Exception get!(X, :x, 1)
    @test_throws Exception pop!(X, :a)
    @test_throws Exception pop!(X, :x, 1)

    # Annotated arrays with mutable properties.
    kwds = (; a=1, b=2, c=3)
    @test @inferred(AnnotatedArrays.keytype(kwds)) === Symbol
    @test @inferred(AnnotatedArrays.keytype(typeof(kwds))) === Symbol
    @test @inferred(AnnotatedArrays.valtype(kwds)) === Int
    @test @inferred(AnnotatedArrays.valtype(typeof(kwds))) === Int
    dict = Dict{Symbol,Int}((k=>v for (k,v) in pairs(kwds))...)
    @test @inferred(AnnotatedArrays.keytype(dict)) === keytype(dict)
    @test @inferred(AnnotatedArrays.keytype(typeof(dict))) === keytype(typeof(dict))
    @test @inferred(AnnotatedArrays.valtype(dict)) === valtype(dict)
    @test @inferred(AnnotatedArrays.valtype(typeof(dict))) === valtype(typeof(dict))
    X = @inferred(AnnotatedArray{T,N}(A, pairs(dict)...))
    @test isa(X, AnnotatedArray{T,N})
    # Accessors.
    @test @inferred(values(X)) === A
    @test @inferred(properties(X)) == dict
    @test keytype(properties(X)) === keytype(dict)
    @test valtype(properties(X)) === Any
    @test length(properties(X)) === length(dict)
    # Other constructors sharing `A` for values.
    Y = @inferred(AnnotatedArray{T}(A, dict))
    @test isa(Y, AnnotatedArray{T,N})
    @test @inferred(values(Y)) === A
    @test @inferred(properties(Y)) === dict
    @test @inferred(AnnotatedArray{T,N}(A, dict)) === Y
    Z = @inferred(AnnotatedArray(A; kwds...))
    @test isa(Z, AnnotatedArray{T,N})
    @test @inferred(values(Z)) === A
    @test @inferred(properties(Z)) == dict
    @test keytype(properties(Z)) === keytype(dict)
    @test valtype(properties(Z)) === Any
    @test length(properties(Z)) === length(dict)
    # Other constructors with a new array for values.
    @test_throws Exception AnnotatedArray(undef, dims, dict) # NOTE missing element type `T`
    Y = @inferred(AnnotatedArray{T}(undef, dims, dict))
    @test isa(Y, AnnotatedArray{T,N})
    @test size(Y) == dims
    @test @inferred(values(Y)) isa Array{T,N}
    @test @inferred(properties(Y)) === dict
    Z = @inferred(AnnotatedArray{T,N}(undef, dims; kwds...))
    @test isa(Z, AnnotatedArray{T,N})
    @test size(Z) == dims
    @test @inferred(values(Z)) isa Array{T,N}
    @test @inferred(properties(Z)) == dict
    # Abstract array API.
    @test @inferred(length(Y)) == length(values(Y)) == length(A)
    @test @inferred(size(Y)) === size(A)
    @test @inferred(axes(Y)) === axes(A)
    @test @inferred(ndims(Y)) === @inferred(ndims(typeof(Y))) === ndims(A)
    @test @inferred(eltype(Y)) === @inferred(eltype(typeof(Y))) === eltype(A)
    @test @inferred(firstindex(Y)) === firstindex(A)
    @test @inferred(lastindex(Y)) === lastindex(A)
    @test @inferred(CartesianIndices(Y)) == CartesianIndices(A)
    @test @inferred(IndexStyle(Y)) === IndexLinear()
    @test @inferred(fill!(Y, zero(T))) === Y
    @test Y != A # only compare values (FIXME?)
    @test @inferred(copyto!(Y, A)) === Y
    @test Y == A # only compare values (FIXME?)
    # Properties.
    @test @inferred(properties(X)) isa Dict{Symbol,Any}
    @test length(properties(X)) == 3
    let names = @inferred(propertynames(X))
        @test length(names) == length(properties(X)) == 3
        @test (:a ∈ names) && (:b ∈ names) && (:c ∈ names)
    end
    @test X.a === kwds.a
    @test X.b === kwds.b
    @test X.c === kwds.c
    @test_throws Exception X.unknown
    # Indexation by array index or by property key.
    i_start, i_stop = firstindex(A), lastindex(A)
    rng = i_start : min(i_start + 3, i_stop)
    I_start, I_stop = first(CartesianIndices(A)), last(CartesianIndices(A))
    I_step = CartesianIndex(ntuple(i -> i == 1 ? 1 : 0, Val(N)))
    @test X[i_start] == A[i_start]
    @test X[i_stop] == A[i_stop]
    @test X[I_start] == A[I_start]
    @test X[I_stop] == A[I_stop]
    @test X[rng] == A[rng]
    @test (X[:a], X[:b], X[:c]) === (kwds.a, kwds.b, kwds.c)
    @test  haskey(X, :a)
    @test !haskey(X, :x)
    @test_throws MethodError haskey(X, i_start) # NOTE `haskey` is not implemented for arrays
    @test getkey(X, :b, undef) === :b
    @test getkey(X, :x, undef) === undef
    @test_throws MethodError getkey(X, i_start, undef) # NOTE `getkey` is not implemented for arrays
    @test get(X, :c, undef) == X.c
    @test get(X, :x, undef) == undef
    @test get(X, firstindex(X), undef) == first(A)
    @test get(X, firstindex(X) - 1, undef) === undef
    @test get(X, lastindex(X), undef) === last(A)
    @test get(X, lastindex(X) + 1, undef) === undef
    @test get(X, I_start, undef) === first(A)
    @test get(X, I_start - I_step, undef) === undef
    @test get(X, I_stop, undef) === last(A)
    @test get(X, I_stop + I_step, undef) === undef
    delete!(X, :a)
    @test !haskey(X, :a)
    @test :a ∉ propertynames(X)
    X.a = 17
    @test haskey(X, :a) && X.a == 17
    @test get!(X, :a, -123) == 17
    @test !haskey(X, :z)
    @test get!(X, :z, -101) == -101
    @test get(X, :z, undef) == -101
    @test pop!(X, :z) == -101
    @test !haskey(X, :z)
    @test_throws Exception pop!(X, :z) # z is not defined
    @test pop!(X, :z, undef) === undef
    @test get(() -> 345, X, :z) == 345
    @test !haskey(X, :z)
    @test get!(() -> 345, X, :z) == 345
    @test haskey(X, :z) && X.z == 345
    @test :z ∈ propertynames(X)

    # Annotated arrays with string keys.
    dict = Dict{String,Any}("units" => "photons", "Δx" => 0.20, "Δy" => 0.15)
    @test @inferred(AnnotatedArrays.keytype(dict)) === keytype(dict)
    @test @inferred(AnnotatedArrays.keytype(typeof(dict))) === keytype(typeof(dict))
    @test @inferred(AnnotatedArrays.valtype(dict)) === valtype(dict)
    @test @inferred(AnnotatedArrays.valtype(typeof(dict))) === valtype(typeof(dict))
    X = @inferred(AnnotatedArray(A, pairs(dict)...))
    @test X isa AnnotatedArray{T,N,Dict{String,Any}}
    @test properties(X) isa Dict{String,Any}
    @test properties(X) !== dict && properties(X) == dict
    names = @inferred propertynames(X)
    @test sort(collect(names)) == sort(collect(keys(dict)))
    @test "units" ∈ names
    @test X["units"] == dict["units"]
    @test "Δx" ∈ names
    @test X["Δx"] == dict["Δx"]
    @test haskey(X, "Δy")

    # Annotated array with no initial properties.
    X = @inferred(AnnotatedArray(A))
    @test properties(X) isa Dict{Symbol,Any}
    @test isempty(propertynames(X))
    X.a = 1
    @test length(properties(X)) == 1
    @test haskey(X, :a)

end

@testset "Pseudo-arrays " begin
#=
    Q = UnfinishedArray(Va)
    @test_throws ErrorException Q[1]
    R = DummyArray(Va)

    copyto!(G, rand(dims...)) # exercise setindex! for linear indices
    copyto!(R, rand(size(Va)...)) # exercise setindex! for Cartesian indices
    sum1 = 0.0
    for val in G
        sum1 += val
    end
    @test samevalues(F, G)
    @test sum(F) == sum(G) ≈ sum1
    sum2 = 0.0
    for val in R
        sum2 += val
    end
    @test samevalues(R, Va)
    @test sum(R) == sum(Va) ≈ sum2
    @test pointer(parent(F)) == pointer(parent(G))
    @test length(F) == length(G) == prod(dims)
    @test eltype(F) == eltype(G) == T
    @test ndims(F) == ndims(G) == N
    @test size(F) == size(G) == dims
    @test ntuple(d -> size(F,d), N) == ntuple(d -> size(G,d), N) == dims
    @test axes(F) == axes(G) == inds
    @test ntuple(d -> axes(F,d), N) == ntuple(d -> axes(G,d), N) == inds
    @test Base.axes1(F) == Base.axes1(G) == inds[1]
    @test IndexStyle(F) == IndexStyle(G) == IndexStyle(parent(F))
    @test_throws ErrorException parent(Q)
    @test IndexStyle(R) == IndexStyle(parent(R)) == IndexCartesian()
    @test pairs(IndexStyle(F), F) == pairs(IndexStyle(G), G) == pairs(IndexStyle(parent(F)), parent(F))
    @test pairs(IndexStyle(R), R) == pairs(IndexStyle(Va), Va)
    @test_throws ErrorException size(F,0)
    @test_throws BoundsError F[0]
    @test keytype(F) === (ndims(F) == 1 ? Int : CartesianIndex{ndims(H)})
    keytype(G) === String
    @test keytype(H) === CartesianIndex{ndims(H)}
    @test valtype(F) === eltype(F)
    @test valtype(G) === eltype(G)
    @test valtype(H) === eltype(H)

    @test length(properties(F)) == length(properties(properties(F))) == 0 && length(properties(G)) == length(properties(properties(G))) == 3
    @test keys(F) == keys(properties(F))
    @test keys(G) == keys(properties(G))
    @test values(G) == values(properties(G))
    @test pairs(G) == pairs(properties(G))
    @test_throws ArgumentError (H["ga"] = π)
    @test haskey(H, "ga") == false
    @test (H[:ga] = π) ≈ π
    @test haskey(H, :ga) == true
    @test keys(H) == keys(properties(H))
    @test getkey(H, :ga, nothing) === :ga
    @test getkey(H, :bu, nothing) === nothing
    @test getkey(H, "ga", nothing) === nothing
    @test get(H, :ga, nothing) ≈ π
    @test get(H, :bu, nothing) === nothing
    @test get(H, "ga", nothing) === nothing
    @test get!(H, :ga, nothing) ≈ π
    @test get!(H, :bu, 55) == 55 && haskey(H, :bu) == true
    @test_throws Exception get!(H, "ga", nothing)
    @test length(properties(F)) == length(properties(G))
    @test keys(F) == keys(G)
    @test length(properties(F)) == length(properties(G)) + 1
    @test pop!(F, "ga") == π
    @test_throws KeyError F["ga"]
    @test_throws KeyError pop!(F, "ga")
    @test pop!(F, "ga", 42) == 42
    @test length(properties(delete!(F, "ga"))) == length(properties(G))
    @test length(properties(delete!(F, "units"))) == length(properties(G)) - 1
    Dat1 = Dict{String,Any}("a" => 1, "b" => sqrt(2), "c" => 4)
    Dat2 = Dict{String,Any}("a" => -1.0, "b" => π, "c" => sqrt(3), "d" => 11)
    for k in keys(F)
        delete!(F, k)
    end
    for k in keys(G)
        delete!(G, k)
    end
=#
end

end # module
