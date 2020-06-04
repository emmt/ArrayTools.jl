module ArrayToolsTests

using Test, Random
using ArrayTools, ArrayTools.AnnotatedArrays, ArrayTools.PseudoArrays
using .AnnotatedArrays:
    propertyname,
    propertynames,
    DynamicallyAnnotatedArray,
    StaticallyAnnotatedArray
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

@testset "Miscellaneous" begin
    # Promotion of element types.
    for (T1,T2,T3) in ((Float32,Float64,Int), (Int16,Int32,Int64))
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
    # Dimensions.
    @test dimensions(()) === ()
    @test dimensions(5) === (5,)
    @test dimensions(Int16(5)) === (5,)
    @test isa(dimensions(Int16(5)), Tuple{Int})
    @test dimensions(A) === size(A)
    @test dimensions(dims) === dims
    @test dimensions(dims...) === dims
    @test dimensions(UInt16(dims[1]), dims[2:end]...) === dims
    @test dimensions((UInt16(dims[1]), dims[2:end]...)) === dims
    @test check_dimensions(()) === true
    @test check_dimensions((1,0,2)) === true
    @test_throws ErrorException check_dimensions((1,0,-1))
    @test_throws ErrorException check_dimensions(-1)
    @test_deprecated checkdimensions((1,2,))
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
    @test_deprecated indices(A) === CartesianIndices(axes(A))
    @test cartesian_indices(A) === CartesianIndices(axes(A))
    @test cartesian_indices(A) === cartesian_indices(size(A))
    @test cartesian_indices(A) === cartesian_indices(axes(A))
    @test cartesian_indices(A) === cartesian_indices(cartesian_indices(A))
    I1 = CartesianIndex(1,2,3)
    I2 = CartesianIndex(5,7,9)
    @test cartesian_indices(I1,I2) ===
        CartesianIndices(([I1[k]:I2[k] for k in 1:length(I1)]...,))
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
    @test all_match(nothing, identity)
    @test all_match(1, length, π)
    @test all_match(1, length, π, 2, "e")
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
    @test_deprecated safe_indices(X) === all_indices(X)
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
    @test_throws ErrorException ArrayTools.throw_non_standard_indexing()
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
    @test_deprecated colons(5) === rubberindex(5)
    #@test_deprecated A[…] === A # FIXME: not possible to deprecate
    for d ∈ 0:12
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

@testset "Storage" begin
    B = flatarray(Float32, A)
    C = flatarray(Float32, Va)
    @test StorageType() === AnyStorage()
    @test StorageType("a") === AnyStorage()
    @test StorageType(A) === FlatStorage()
    @test StorageType(Vf) === FlatStorage()
    @test StorageType(Va) === AnyStorage()
    @test isflatarray() == false
    @test isflatarray("a") == false
    @test isflatarray(S) == false
    @test isflatarray(A) == true
    @test isflatarray(Vf) == true
    @test isflatarray(Va) == false
    @test isflatarray(A,B,C) == true
    @test isflatarray(A,B,C,Vf) == true
    @test isflatarray(A,B,C,Va) == false
    @test flatarray(A) === A
    @test flatarray(eltype(A), A) === A
    @test pointer(A) != pointer(flatarray(Va))
    @test samevalues(S, flatarray(S))
    @test eltype(S) === eltype(flatarray(S))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(Va, C) ≤ atol
    Da = flatarray(Va)
    @test isflatarray(Da) == true
    @test Da == Va
    @test Da !== Va
    Df = flatarray(Vf)
    @test isflatarray(Df) == true
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
        @test isflatarray(x) == b
    end
    @test isflatarray(A, view(A, :, :, 2, 3), view(A, :, :, 2:2, 3)) == true
    @test isflatarray(A, view(A, :, :, 2:2, :), view(A, :, :, 2:2, 3)) == false
end

@testset "Indexing" begin
    for Q in (A,Va,S,101,Colon())
        @test has_standard_indexing(Q) == !Base.has_offset_axes(Q)
    end
    B = fastarray(Float32, A)
    C = fastarray(Float32, Va)
    @test has_standard_indexing(A,Va) == (has_standard_indexing(A) &&
                                         has_standard_indexing(Va))
    @test IndexingTrait(A) === FastIndexing()
    @test IndexingTrait(Va) === AnyIndexing()
    @test IndexingTrait("a") === AnyIndexing()
    @test isfastarray() == false
    @test isfastarray(S) == true
    @test samevalues(Va, fastarray(Va))
    @test samevalues(A, fastarray(A))
    @test isfastarray(A) == true
    @test isfastarray(A,B,C) == true
    @test isfastarray(Va) == false
    @test isfastarray(A,B,C,Va) == false
    @test isfastarray(fastarray(Va)) == true
    @test pointer(A) == pointer(fastarray(A))
    @test pointer(A) != pointer(fastarray(Va))
    @test samevalues(Va, fastarray(Va))
    @test samevalues(A, fastarray(A))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(Va, C) ≤ atol
end

@testset "Broadcasting" begin
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

@testset "Custom arrays" begin
    inds = map(n -> Base.OneTo(n), dims)
    T = Float32
    N = length(dims)
    D1 = Dict("units" => "photons", "Δx" => 0.20, "Δy" => 0.15)
    @test isa(D1, Dict{String,Any})
    K2 = (x=true, y=1.8, units="µm")
    D2 = Dict(pairs(K2)...)
    @test isa(D2, Dict{Symbol,Any})
    D3 = Dict(:x => 1, :y => 2, :z => 3)
    @test isa(D3, Dict{Symbol,Int})
    D4 = Dict(:x => 1, Float64 => 3.1, "string" => "hello")
    @test isa(D4, Dict{Any,Any})
    D5 = Dict(1.0 => 1, 2.0 => :2, 3.0 => "3")
    @test isa(D5, Dict{Float64,Any})
    G = AnnotatedArray(zeros(T, dims), pairs(D1)...)
    F = AnnotatedArray{T}(parent(G), Dict{String,Any}())
    H = AnnotatedArray{T,N}(undef, dims, Dict{Symbol,Float32}())
    @test (AnnotatedArrays.keytype(D1) ===
           AnnotatedArrays.keytype(typeof(D1)) ===
           keytype(D1) === keytype(typeof(D1)) === String)
    @test (AnnotatedArrays.valtype(D1) ===
           AnnotatedArrays.valtype(typeof(D1)) ===
           valtype(D1) === valtype(typeof(D1)) === Any)
    @test (AnnotatedArrays.keytype(D2) ===
           AnnotatedArrays.keytype(typeof(D2)) ===
           keytype(D2) === keytype(typeof(D2)) === Symbol)
    @test (AnnotatedArrays.valtype(D2) ===
           AnnotatedArrays.valtype(typeof(D2)) ===
           valtype(D2) === valtype(typeof(D2)) === Any)
    @test (AnnotatedArrays.keytype(K2) ===
           AnnotatedArrays.keytype(typeof(K2)) === Symbol)
    @test (AnnotatedArrays.valtype(K2) ===
           AnnotatedArrays.valtype(typeof(K2)) === Any)
    @test (AnnotatedArrays.keytype(D3) ===
           AnnotatedArrays.keytype(typeof(D3)) ===
           keytype(D3) === keytype(typeof(D3)) === Symbol)
    @test (AnnotatedArrays.valtype(D3) ===
           AnnotatedArrays.valtype(typeof(D3)) ===
           valtype(D3) === valtype(typeof(D3)) === Int)
    @test (AnnotatedArrays.keytype(D4) ===
           AnnotatedArrays.keytype(typeof(D4)) ===
           keytype(D4) === keytype(typeof(D4)) === Any)
    @test (AnnotatedArrays.valtype(D4) ===
           AnnotatedArrays.valtype(typeof(D4)) ===
           valtype(D4) === valtype(typeof(D4)) === Any)
    @test (AnnotatedArrays.keytype(D5) ===
           AnnotatedArrays.keytype(typeof(D5)) ===
           keytype(D5) === keytype(typeof(D5)) === Float64)
    @test (AnnotatedArrays.valtype(D5) ===
           AnnotatedArrays.valtype(typeof(D5)) ===
           valtype(D5) === valtype(typeof(D5)) === Any)

    # Forbidden key types.
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{Any,Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{Int32,Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{CartesianIndex,Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{UnitRange{Int},Any}())
    @test_throws ErrorException AnnotatedArray{T}(undef, dims, Dict{Colon,Any}())

    # Common errors.
    @test_throws ErrorException AnnotatedArray(undef, dims)
    @test_throws ErrorException AnnotatedArray(undef, dims...)

    # Annotated arrays with immutable properties.
    let X = AnnotatedArray(A, (a = 1,)),
        Y = AnnotatedArray{T}(A, (a = "1", b = "2"))
        @test isa(X, StaticallyAnnotatedArray)
        @test (AnnotatedArrays.keytype(X) === keytype(X) ===
               AnnotatedArrays.keytype(typeof(X)) === keytype(typeof(X)) === Symbol)
        @test (AnnotatedArrays.valtype(X) === valtype(X) ===
               AnnotatedArrays.valtype(typeof(X)) === valtype(typeof(X)) === Int)
        @test nkeys(X) == 1
        @test propertyname(typeof(X), :foo) === :foo
        @test isa(Y, StaticallyAnnotatedArray)
        @test (AnnotatedArrays.keytype(Y) === keytype(Y) ===
               AnnotatedArrays.keytype(typeof(Y)) === keytype(typeof(Y)) === Symbol)
        @test (AnnotatedArrays.valtype(Y) === valtype(Y) ===
               AnnotatedArrays.valtype(typeof(Y)) === valtype(typeof(Y)) === String)
        @test nkeys(Y) == 2
        @test_throws ErrorException X.field = 40
        @test_throws ErrorException X[:field] = 40
        @test :a ∈ propertynames(X)
        @test :b ∈ propertynames(Y)
        @test haskey(X, :a) == true
        @test haskey(X, :b) == false
        @test haskey(X, "b") == false
        @test getkey(X, :a, :x) == :a
        @test getkey(X, :b, :x) == :x
        @test getkey(X, "b", :x) == :x
        @test get(X, :a, :x) == X.a
        @test get(X, :b, :x) == :x
        @test get(X, "b", :x) == :x
        @test X[:a] == X.a
        @test_throws ErrorException X.unknown
        @test_throws ErrorException delete!(X, :a)
        @test_throws ErrorException get!(X, :a, 1)
        @test_throws ErrorException pop!(X, :a)
        @test_throws ErrorException pop!(X, :a, 1)
        @test_throws ErrorException merge!(X, Y)
        @test_throws ErrorException merge!(+, X, Y)
    end
    let X = AnnotatedArray{T}(undef, dims, (a = 1.0, b = 2.0, c = 3.0))
        @test isa(X, StaticallyAnnotatedArray)
        @test eltype(X) === T
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Float64
        @test nkeys(X) == 3
    end
    let X = AnnotatedArray{T,N}(A, (a = 1, b = 2, c = 3, d = 4))
        @test isa(X, StaticallyAnnotatedArray)
        @test eltype(X) === T
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Int
        @test nkeys(X) == 4
    end
    let X = AnnotatedArray{T,N}(undef, dims,
                                (a = 1, b = 2, c = 3, d = 4, e = 5))
        @test isa(X, StaticallyAnnotatedArray)
        @test eltype(X) === T
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Int
        @test nkeys(X) == 5
    end

    # Annotated arrays with mutable properties.
    let X = AnnotatedArray{T,N}(A)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T,N}(A, D3)
        @test isa(X, DynamicallyAnnotatedArray)
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Int
    end
    let X = AnnotatedArray{T,N}(A, pairs(D3)...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    @test_throws ErrorException AnnotatedArray{T,N}(A, D4)
    @test_throws ErrorException AnnotatedArray{T,N}(A, pairs(D4))
    @test_throws ErrorException AnnotatedArray{T,N}(A, pairs(D4)...)
    let X = AnnotatedArray{T,N}(undef, dims)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
        X.foo = 42
        X.bar = π
        @test X.foo === X[:foo] === 42
        @test X.bar === X[:bar] === π
    end
    let X = AnnotatedArray{T,N}(undef, dims, D2)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T,N}(undef, dims, pairs(D2)...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T,N}(A; K2...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T,N}(undef, dims; K2...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(A)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(A, D1)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === String
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
        @test propertyname(typeof(X), :somekey) == "somekey"
    end
    let X = AnnotatedArray{T}(A, pairs(D1)...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === String
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(undef, dims)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(undef, dims, D2)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(undef, dims, pairs(D2)...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(undef, dims; K2...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray(A)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray(A, D1)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === String
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray(A, pairs(D2)...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray(A; K2...)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test ndims(X) == N
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray(A, π => "π", sqrt(2) => "√2", 1.3 => :bar)
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === eltype(A)
        @test axes(X) === axes(A)
        @test nkeys(X) == 3
        @test AnnotatedArrays.keytype(X) === keytype(X) === Real
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
        @test_throws ErrorException propertyname(typeof(X), :foo)
    end
    let X = AnnotatedArray{T}(undef, dims, π=>"π", sqrt(2)=>"√2")
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test size(X) === dims
        @test nkeys(X) == 2
        @test AnnotatedArrays.keytype(X) === keytype(X) === Real
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
    end
    let X = AnnotatedArray{T}(undef, dims...; a="π", b="√2")
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
        @test propertyname(typeof(X), :foo) === :foo
    end
    let X = AnnotatedArray{T,N}(undef, dims...; a="π", b="√2")
        @test isa(X, DynamicallyAnnotatedArray)
        @test eltype(X) === T
        @test size(X) === dims
        @test AnnotatedArrays.keytype(X) === keytype(X) === Symbol
        @test AnnotatedArrays.valtype(X) === valtype(X) === Any
        @test X.a == "π"
        @test X.b == "√2"
        @test :a ∈ propertynames(X)
    end

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
    @test Base.elsize(F) == Base.elsize(G) == Base.elsize(parent(F))
    @test Base.elsize(R) == Base.elsize(typeof(R)) == Base.elsize(Va)
    @test sizeof(F) == sizeof(G) == sizeof(parent(F))
    @test IndexStyle(F) == IndexStyle(G) == IndexStyle(parent(F))
    @test_throws ErrorException parent(Q)
    @test IndexStyle(R) == IndexStyle(parent(R)) == IndexCartesian()
    @test pairs(IndexStyle(F), F) == pairs(IndexStyle(G), G) == pairs(IndexStyle(parent(F)), parent(F))
    @test pairs(IndexStyle(R), R) == pairs(IndexStyle(Va), Va)
    @test_throws ErrorException size(F,0)
    @test_throws BoundsError F[0]
    @test keytype(F) === keytype(G) === String
    @test keytype(H) === Symbol
    @test valtype(F) === valtype(G) === Any
    @test valtype(H) === Float32

    @test nkeys(F) == nkeys(properties(F)) == 0 && nkeys(G) == nkeys(properties(G)) == 3
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
    merge!(F, G)
    @test nkeys(F) == nkeys(G)
    @test keys(F) == keys(G)
    merge!(F, G, Dict("ga" => π))
    @test nkeys(F) == nkeys(G) + 1
    @test pop!(F, "ga") == π
    @test_throws KeyError F["ga"]
    @test_throws KeyError pop!(F, "ga")
    @test pop!(F, "ga", 42) == 42
    @test nkeys(delete!(F, "ga")) == nkeys(G)
    @test nkeys(delete!(F, "units")) == nkeys(G) - 1
    @test merge(G, F) == merge(properties(G), properties(F))
    @test merge(D2, F, D1) == merge(D2, properties(F), D1)
    Dat1 = Dict{String,Any}("a" => 1, "b" => sqrt(2), "c" => 4)
    Dat2 = Dict{String,Any}("a" => -1.0, "b" => π, "c" => sqrt(3), "d" => 11)
    for k in keys(F)
        delete!(F, k)
    end
    for k in keys(G)
        delete!(G, k)
    end
    merge!(F, Dat1)
    merge!(G, Dat2)
    @test merge(+, F, Dat2) == merge(+, Dat1, Dat2)
    @test merge(-, Dat1, G) == merge(-, Dat1, Dat2)
    @test merge(*, F, G) == merge(*, Dat1, Dat2)
    @test properties(merge!(+, F, Dat2)) == merge(+, Dat1, Dat2)
    @test merge!(+, copy(Dat1), G) == properties(F)
end

end # module
