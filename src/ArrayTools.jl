module ArrayTools

export
    ..,
    ArraySize,
    RubberIndex,
    all_match,
    allof,
    anyof,
    axis_limits,
    bcastcopy,
    bcastlazy,
    bcastsize,
    cartesian_indices,
    check_size,
    colons,
    common_indices,
    has_standard_indexing,
    noneof,
    promote_eltype,
    reversemap,
    all_indices,
    same_axes,
    same_size,
    same_standard_size,
    split_interval,
    standard_size,
    strictmap!,
    to_int,
    to_size,
    # storage trait
    StorageType,
    AnyStorage,
    FlatStorage,
    to_flat_array,
    is_flat_array,
    # Fast arrays and indexing trait:
    IndexingType,
    FastIndexing,
    AnyIndexing,
    to_fast_array,
    is_fast_array

using Base: OneTo, axes1, tail
import Base: dotview, getindex, setindex!, to_indices

@deprecate rubberindex colons
@deprecate indices cartesian_indices
@deprecate cartesianindices cartesian_indices
@deprecate safe_indices all_indices
@deprecate flatmatrix to_flat_array
@deprecate flatvector to_flat_array
@deprecate fastmatrix to_fast_array
@deprecate fastvector to_fast_array
@deprecate isfastarray is_fast_array
@deprecate isflatarray is_flat_array
@deprecate fastarray to_fast_array
@deprecate flatarray to_flat_array
@deprecate dimensions(A::AbstractArray) standard_size(A)
@deprecate dimensions(dim::Integer) to_size(dims)
@deprecate dimensions(dims::Integer...) to_size(dims)
@deprecate dimensions(dims::Tuple{Vararg{Integer}}) to_size(dims)
@deprecate bcastdim(a::Integer, b::Integer) bcastsize(a, b)
@deprecate bcastdims bcastsize
@deprecate Dimensions ArraySize
@deprecate checkdimensions check_size
@deprecate check_dimensions check_size
@deprecate same_dimensions same_standard_size

include("traits.jl")
include("utils.jl")
include("indexing.jl")
include("broadcasting.jl")
include("rubberindex.jl")

include("PseudoArrays.jl")
using .PseudoArrays

include("AnnotatedArrays.jl")
using .AnnotatedArrays

end # module
