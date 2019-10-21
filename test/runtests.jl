module ArrayToolsTests
using Test, ArrayTools

function same(A::AbstractArray, B::AbstractArray)
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

dims = (3, 4, 5)
A = rand(Float64, dims)
V = view(A, :, 2:3, :)
R = 1:2:100
atol = 1e-6

@testset "Miscellaneous" begin
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
    # Colons.
    for d ∈ 0:12
        tup = ntuple(x -> Colon(), d)
        @test colons(d) === tup
        @test colons(Val(d)) === tup
    end
    for k ∈ 1:dims[end]
        @test same(slice(A, k), A[:,:,k])
    end
    # Other stuff.
    @test reversemap(x -> x^2, dims) === reverse(map(x -> x^2, dims))
end

@testset "Storage" begin
    B = flatarray(Float32, A)
    C = flatarray(Float32, V)
    @test isflatarray() == false
    @test isflatarray(R) == false
    @test isflatarray(A) == true
    @test isflatarray(V) == false
    @test isflatarray(A,B,C) == true
    @test isflatarray(A,B,C,V) == false
    @test isflatarray(flatarray(V)) == true
    @test pointer(A) == pointer(flatarray(A))
    @test pointer(A) == pointer(flatarray(eltype(A), A))
    @test pointer(A) != pointer(flatarray(V))
    @test same(R, flatarray(R))
    @test eltype(R) === eltype(flatarray(R))
    @test same(V, flatarray(V))
    @test same(A, flatarray(A))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(V, C) ≤ atol
    for n in 1:5
        K = rand(Float64, ntuple(x -> 3, n))
        L = view(K, 2:3, colons(n-1)...)
        @test StorageType(K) == FlatStorage()
        @test StorageType(L) == (n == 1 ? FlatStorage() : AnyStorage())
    end
end

@testset "Indexing" begin
    for Q in (A,V,R)
        @test has_standard_indexing(Q) == !Base.has_offset_axes(Q)
    end
    B = fastarray(Float32, A)
    C = fastarray(Float32, V)
    @test has_standard_indexing(A,V) == (has_standard_indexing(A) &&
                                         has_standard_indexing(V))
    @test IndexingTrait(A) == FastIndexing()
    @test IndexingTrait(V) == AnyIndexing()
    @test isfastarray(R) == true
    @test same(V, fastarray(V))
    @test same(A, fastarray(A))
    @test isfastarray(A) == true
    @test isfastarray(A,B,C) == true
    @test isfastarray(V) == false
    @test isfastarray(A,B,C,V) == false
    @test isfastarray(fastarray(V)) == true
    @test pointer(A) == pointer(fastarray(A))
    @test pointer(A) != pointer(fastarray(V))
    @test same(V, fastarray(V))
    @test same(A, fastarray(A))
    @test maxabsdif(A, B) ≤ atol
    @test maxabsdif(V, C) ≤ atol
end

end # module
