module ArrayToolsTests

using Test, Random
using ArrayTools, ArrayTools.AnnotatedArrays, ArrayTools.PseudoArrays

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

T = Float32
dims = (3, 4, 5, 6)
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
    # Tests for `safe_indices`.
    #
    B = rand(T, dims[1:2])
    C = rand(T, dims[1:end-2]..., dims[end], dims[end-1])
    X = rand(T, 6)
    @test safe_indices(X) === eachindex(X)
    @test safe_indices(A) === eachindex(A)
    @test safe_indices(Va) === eachindex(Va)
    @test_throws DimensionMismatch safe_indices(A,Va)
    @test eachindex(A,C) === eachindex(A)
    @test_throws DimensionMismatch safe_indices(A,C)
    @test safe_indices(A, rand(T, dims), rand(T, dims)) === eachindex(A)
    Y = rand(T, dims[1], dims[2]+2, dims[3:end]...)
    Z = view(Y, :, 2:dims[2]+1, colons(length(dims)-2)...)
    @test IndexStyle(Z) === IndexCartesian()
    @test safe_indices(A, rand(T, dims), Z) === eachindex(IndexCartesian(), A)
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
    for d ∈ 0:12
        tup = ntuple(x -> Colon(), d)
        @test colons(d) === tup
        @test colons(Val(d)) === tup
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
    # The `end` keyword can only appear **before** the rubber index.
    @test_broken A[2,…,3:end] == A[2,:,:,3:end]
    A[1,…] .= 3
    @test all(isequal(3), A[1,:,:,:])
    A[1,…] = A[2,…]
    @test A[1,:,:,:] == A[2,:,:,:]
    A[1,…,2:3] .= A[2,…,3:2:5]
    @test A[1,:,:,2:3] == A[2,:,:,3:2:5]
    # Tests with a 9 dimensional array.
    siz = (3, 4, 2, 3, 4, 2, 3, 4, 2)
    n = length(siz)
    B = generate(Int, siz)
    C = copy(B)
    inds = (2:3, 4, Colon(), CartesianIndex(2),
            CartesianIndex(3), 1:2, 3, 2:3, CartesianIndex(1))
    for k ∈ 0:n
        @test ArrayTools.numberofindices(inds[1:k]...) == k
    end
    for k in 0:n, l in 0:n-k
        J1 = inds[1:k]         # leading indices
        J2 = inds[n-l+1:n]     # trailing indices
        J = (J1..., …, J2...,) # with a rubber index
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
    @test_throws DimensionMismatch bcastdim(3, 4)
    for (a,b) in (((), ()),
                  ((Int16(5), Int32(6), Int64(7)), (5,6,7)),)
        @test bcastdims(a) === b
        @test bcastdims((), a) === b
        @test bcastdims(a, ()) === b
        @test bcastdims(a, (a..., 3,)) === (b..., 3,)
        @test bcastdims(a, (a..., 3,), (a..., 1, 2)) === (b..., 3, 2)
        @test bcastdims((a..., 3,), a) === (b..., 3,)
        @test bcastdims((a..., 1,), (a..., 3,)) === (b..., 3,)
        @test bcastdims((a..., 3,), (a..., 1,)) === (b..., 3,)
    end
    # Check `bcastcopy`.
    T = Int32
    dims1 = (1,2,3)
    dims2 = (3,2,1,4)
    dims3 = bcastdims(dims1, dims2)
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
    B4 = bcastlazy(Float64, 1, dims...)
    @test eltype(B4) === Float64
    @test size(B4) === dims
    @test B4 == ones(Float64, dims)
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
    @test_throws ErrorException delete!(A11, :a)
    @test_throws ErrorException get!(A11, :a, 1)
    @test_throws ErrorException pop!(A11, :a)
    @test_throws ErrorException pop!(A11, :a, 1)
    @test_throws ErrorException merge!(A11, A12)
    @test_throws ErrorException merge!(+, A11, A12)

    A31 = AnnotatedArray{T,N}(Array{T,N}(undef, dims))
    A32 = AnnotatedArray{T,N}(Array{T,N}(undef, dims), D3)
    A33 = AnnotatedArray{T,N}(Array{T,N}(undef, dims), pairs(D3)...)
    @test_throws ErrorException AnnotatedArray{T,N}(Array{T,N}(undef, dims),
                                                    D4)
    @test_throws ErrorException AnnotatedArray{T,N}(Array{T,N}(undef, dims),
                                                    pairs(D4))
    @test_throws ErrorException AnnotatedArray{T,N}(Array{T,N}(undef, dims),
                                                    pairs(D4)...)
    A34 = AnnotatedArray{T,N}(undef, dims)
    A35 = AnnotatedArray{T,N}(undef, dims, D2)
    A36 = AnnotatedArray{T,N}(undef, dims, pairs(D2)...)
    A37 = AnnotatedArray{T,N}(Array{T,N}(undef, dims); K2...)
    A38 = AnnotatedArray{T,N}(undef, dims; K2...)

    A41 = AnnotatedArray{T}(Array{T,N}(undef, dims))
    A42 = AnnotatedArray{T}(Array{T,N}(undef, dims), D1)
    A43 = AnnotatedArray{T}(Array{T,N}(undef, dims), pairs(D1)...)
    A44 = AnnotatedArray{T}(undef, dims)
    A45 = AnnotatedArray{T}(undef, dims, D2)
    A46 = AnnotatedArray{T}(undef, dims, pairs(D2)...)
    A47 = AnnotatedArray{T}(undef, dims; K2...)

    A51 = AnnotatedArray(Array{T,N}(undef, dims))
    A52 = AnnotatedArray(Array{T,N}(undef, dims), D1)
    A53 = AnnotatedArray(Array{T,N}(undef, dims), pairs(D2)...)
    A54 = AnnotatedArray(Array{T,N}(undef, dims); K2...)

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
