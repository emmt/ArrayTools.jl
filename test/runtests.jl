module ArrayToolsTests
using Test, Random, ArrayTools, ArrayTools.AnnotatedArrays, ArrayTools.PseudoArrays

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
    A[rubberindex(Val(N-1))..., i]

dims = (3, 4, 5, 6)
A = rand(Float64, dims)
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
    @test dimensions(()) == ()
    @test dimensions(5) == (5,)
    @test dimensions(Int16(5)) == (5,)
    @test isa(dimensions(Int16(5)), Tuple{Int})
    @test dimensions(A) == size(A)
    @test dimensions(dims) == dims
    @test dimensions(dims...) == dims
    @test dimensions(UInt16(dims[1]), dims[2:end]...) == dims
    @test dimensions((UInt16(dims[1]), dims[2:end]...)) == dims
    @test checkdimensions(()) === true
    @test checkdimensions((1,0,2)) === true
    @test_throws ErrorException checkdimensions((1,0,-1))
    @test_throws ErrorException checkdimensions(-1)
    #
    # Tests for `allof`, `anyof` and `noneof`.
    #
    W = (true, true, true)
    X = [true, false]
    Y = (false, false)
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
    @test allof(collect(W)) == allof(W)
    @test anyof(collect(W)) == anyof(W)
    @test noneof(collect(W)) == noneof(W)
    @test allof(X) == false
    @test anyof(X) == true
    @test noneof(X) == false
    @test allof(Y) == false
    @test anyof(Y) == false
    @test noneof(Y) == true
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
    # Tests for `cartesianindices`.
    #
    @test_deprecated indices(A) === CartesianIndices(axes(A))
    @test cartesianindices(A) === CartesianIndices(axes(A))
    @test cartesianindices(A) === cartesianindices(size(A))
    @test cartesianindices(A) === cartesianindices(axes(A))
    @test cartesianindices(A) === cartesianindices(cartesianindices(A))
    #
    # Tests for `safeindices`.
    #
    B = rand(dims[1:2])
    C = rand(dims[1:end-2]..., dims[end], dims[end-1])
    X = rand(6)
    @test safeindices(X) === eachindex(X)
    @test safeindices(A) === eachindex(A)
    @test safeindices(V) === eachindex(V)
    @test_throws DimensionMismatch safeindices(A,V)
    @test eachindex(A,C) === eachindex(A)
    @test_throws DimensionMismatch safeindices(A,C)
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
    I2 = CartesianIndex(3,4)
    I3 = CartesianIndex(3,4,5)
    @test_deprecated colons(5) == rubberindex(5)
    for d ∈ 0:12
        tup = ntuple(x -> Colon(), d)
        @test rubberindex(d) === tup
        @test rubberindex(Val(d)) === tup
    end
    for k ∈ 1:dims[end]
        @test samevalues(slice(A, k), A[:,:,:,k])
    end
    @test A[…]           == A[:,:,:,:]
    @test A[…,3]         == A[:,:,:,3]
    @test A[2,…]         == A[2,:,:,:]
    @test A[…,2:4,5]     == A[:,:,2:4,5]
    @test A[2:3,…,1,2:4] == A[2:3,:,1,2:4]
    @test A[:,2,…,2:4]   == A[:,2,:,2:4]
    @test A[:,2:3,…]     == A[:,2:3,:,:]
    @test A[…,2:3,:]     == A[:,:,2:3,:]
    @test A[I1,…]        == A[I1,:,:,:]
    @test A[I2,…]        == A[I2,:,:]
    @test A[…,I1]        == A[:,:,:,I1]
    @test A[…,I2]        == A[:,:,I2]
    @test A[I1,…,I2]     == A[I1,:,I2]
    @test A[I2,…,I1]     == A[I2,:,I1]
    @test A[I1,…,I3]     == A[I1,I3]
    @test A[1,2:end,…]   == A[1,2:end,:,:]
    @test A[1,2:end-1,…] == A[1,2:end-1,:,:]
    @test_broken A[2,…,3:end] == A[2,:,:,3:end]
    A[1,…] .= 3
    @test all(isequal(3), A[1,:,:,:])
    A[1,…] = A[2,…]
    @test A[1,:,:,:] == A[2,:,:,:]
    A[1,…,2:3] .= A[2,…,3:2:5]
    @test A[1,:,:,2:3] == A[2,:,:,3:2:5]
end

@testset "Storage" begin
    B = flatarray(Float32, A)
    C = flatarray(Float32, V)
    @test StorageType() == AnyStorage()
    @test StorageType("a") == AnyStorage()
    @test StorageType(A) == FlatStorage()
    @test isflatarray() == false
    @test isflatarray("a") == false
    @test isflatarray(S) == false
    @test isflatarray(A) == true
    @test isflatarray(V) == false
    @test isflatarray(A,B,C) == true
    @test isflatarray(A,B,C,V) == false
    @test isflatarray(flatarray(V)) == true
    @test pointer(A) == pointer(flatarray(A))
    @test pointer(A) == pointer(flatarray(eltype(A), A))
    @test pointer(A) != pointer(flatarray(V))
    @test samevalues(S, flatarray(S))
    @test eltype(S) === eltype(flatarray(S))
    @test samevalues(V, flatarray(V))
    @test samevalues(A, flatarray(A))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(V, C) ≤ atol
    for n in 1:5
        K = rand(Float64, ntuple(x -> 3, n))
        L = view(K, 2:3, rubberindex(n-1)...)
        @test StorageType(K) == FlatStorage()
        @test StorageType(L) == (n == 1 ? FlatStorage() : AnyStorage())
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
    for Q in (A,V,S,101,Colon())
        @test has_standard_indexing(Q) == !Base.has_offset_axes(Q)
    end
    B = fastarray(Float32, A)
    C = fastarray(Float32, V)
    @test has_standard_indexing(A,V) == (has_standard_indexing(A) &&
                                         has_standard_indexing(V))
    @test IndexingTrait(A) == FastIndexing()
    @test IndexingTrait(V) == AnyIndexing()
    @test IndexingTrait("a") == AnyIndexing()
    @test isfastarray() == false
    @test isfastarray(S) == true
    @test samevalues(V, fastarray(V))
    @test samevalues(A, fastarray(A))
    @test isfastarray(A) == true
    @test isfastarray(A,B,C) == true
    @test isfastarray(V) == false
    @test isfastarray(A,B,C,V) == false
    @test isfastarray(fastarray(V)) == true
    @test pointer(A) == pointer(fastarray(A))
    @test pointer(A) != pointer(fastarray(V))
    @test samevalues(V, fastarray(V))
    @test samevalues(A, fastarray(A))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(V, C) ≤ atol
end

@testset "Broadcasting" begin
    A1 = bcastlazy(A, size(A))
    A2 = bcastcopy(A, size(A))
    @test pointer(A1) == pointer(A) && samevalues(A, A1)
    @test pointer(A2) != pointer(A) && samevalues(A, A2)
    A3 = bcastlazy(A, eltype(A), size(A))
    A4 = bcastlazy(A, Float32, size(A))
    @test pointer(A3) == pointer(A) && samevalues(A, A1)
    @test size(A4) == size(A) && maxabsdif(A, A4) ≤ atol
    bdims = (dims[1], 1, dims[3:end]...)
    B = rand(Float64, bdims)
    B1 = bcastlazy(B, dims)
    B2 = bcastcopy(B, dims)
    for i in 1:dims[2]
        @test samevalues(B[:,1,:,:], B1[:,i,:,:])
        @test samevalues(B[:,1,:,:], B2[:,i,:,:])
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
    D2 = Dict(:x => true, :y => 1.8, :units => "µm")
    D3 = Dict(:x => 1, :y => 2, :z => 3)
    G = AnnotatedArray(zeros(T, dims), pairs(D1)...)
    F = AnnotatedArray{T}(parent(G), Dict{String,Any}())
    H = AnnotatedArray{T,N}(undef, dims, Dict{Symbol,Float32}())

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
    A11 = AnnotatedArray(zeros(T, dims), (a = 1,))
    A12 = AnnotatedArray{T}(zeros(T, dims), (a = "1", b = "2"))
    A13 = AnnotatedArray{T}(undef, dims, (a = 1.0, b = 2.0, c = 3.0))
    A14 = AnnotatedArray{T,N}(zeros(T, dims), (a = 1, b = 2, c = 3, d = 4))
    A15 = AnnotatedArray{T,N}(undef, dims, (a = 1, b = 2, c = 3, d = 4, e = 5))
    @test keytype(A11) === keytype(A12) === keytype(A13) === Symbol
    @test valtype(A11) === Int
    @test valtype(A12) === String
    @test valtype(A13) === Float64
    @test nkeys(A11) == 1 && nkeys(A12) == 2 && nkeys(A13) == 3
    @test_throws ErrorException A11.x = 40
    @test_throws ErrorException A11[:x] = 40
    @test haskey(A11, :a) == true
    @test haskey(A11, :b) == false
    @test haskey(A11, "b") == false
    @test getkey(A11, :a, :x) == :a
    @test getkey(A11, :b, :x) == :x
    @test getkey(A11, "b", :x) == :x
    @test get(A11, :a, :x) == A11.a
    @test get(A11, :b, :x) == :x
    @test get(A11, "b", :x) == :x
    @test A11[:a] == A11.a
    @test_throws ErrorException A11.x
    @test_throws ErrorException get!(A11, :a, 1)
    @test_throws ErrorException pop!(A11, :a)
    @test_throws ErrorException pop!(A11, :a, 1)
    @test_throws ErrorException merge!(A11, A12)
    @test_throws ErrorException merge!(+, A11, A12)


    A31 = AnnotatedArray{T,N}(Array{T,N}(undef, dims))
    A32 = AnnotatedArray{T,N}(Array{T,N}(undef, dims), D3)
    A33 = AnnotatedArray{T,N}(Array{T,N}(undef, dims), pairs(D3)...)
    A34 = AnnotatedArray{T,N}(undef, dims)
    A35 = AnnotatedArray{T,N}(undef, dims, D2)
    A36 = AnnotatedArray{T,N}(undef, dims, pairs(D2)...)

    A41 = AnnotatedArray{T}(Array{T,N}(undef, dims))
    A42 = AnnotatedArray{T}(Array{T,N}(undef, dims), D1)
    A43 = AnnotatedArray{T}(Array{T,N}(undef, dims), pairs(D1)...)
    A44 = AnnotatedArray{T}(undef, dims)
    A45 = AnnotatedArray{T}(undef, dims, D2)
    A46 = AnnotatedArray{T}(undef, dims, pairs(D2)...)

    A51 = AnnotatedArray(Array{T,N}(undef, dims))
    A52 = AnnotatedArray(Array{T,N}(undef, dims), D1)
    A53 = AnnotatedArray(Array{T,N}(undef, dims), pairs(D2)...)

    Q = UnfinishedArray(V)
    @test_throws ErrorException Q[1]
    R = DummyArray(V)

    copyto!(G, rand(dims...)) # exercise setindex! for linear indices
    copyto!(R, rand(size(V)...)) # exercise setindex! for Cartesian indices
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
    @test samevalues(R, V)
    @test sum(R) == sum(V) ≈ sum2
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
    @test Base.elsize(R) == Base.elsize(V)
    @test sizeof(F) == sizeof(G) == sizeof(parent(F))
    @test IndexStyle(F) == IndexStyle(G) == IndexStyle(parent(F))
    @test_throws ErrorException parent(Q)
    @test IndexStyle(R) == IndexStyle(parent(R)) == IndexCartesian()
    @test pairs(IndexStyle(F), F) == pairs(IndexStyle(G), G) == pairs(IndexStyle(parent(F)), parent(F))
    @test pairs(IndexStyle(R), R) == pairs(IndexStyle(V), V)
    @test_throws ErrorException size(F,0)
    @test_throws BoundsError F[0]
    @test keytype(F) === keytype(G) === String
    @test keytype(H) === Symbol
    @test valtype(F) === valtype(G) === Any
    @test valtype(H) === Float32
    @test keytype(A31) === Symbol && valtype(A31) == Any
    @test keytype(A32) === Symbol && valtype(A32) == Int
    @test keytype(A33) === Symbol && valtype(A33) == Any

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
