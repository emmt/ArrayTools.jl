module ArrayTools

export
    ..,
    Dimensions,
    RubberIndex,
    all_match,
    allof,
    anyof,
    axis_limits,
    bcastcopy,
    bcastlazy,
    bcastsize,
    cartesian_indices,
    check_dimensions,
    colons,
    common_indices,
    dimensions,
    has_standard_indexing,
    noneof,
    promote_eltype,
    reversemap,
    all_indices,
    same_dimensions,
    same_axes,
    split_interval,
    strictmap!,
    to_int,
    # storage trait
    StorageType,
    AnyStorage,
    FlatStorage,
    flatarray,
    isflatarray,
    # Fast arrays and indexing trait:
    IndexingTrait,
    FastIndexing,
    AnyIndexing,
    fastarray,
    isfastarray

using Base: OneTo, axes1, tail
import Base: dotview, getindex, setindex!, to_indices

@deprecate rubberindex colons
@deprecate indices cartesian_indices
@deprecate cartesianindices cartesian_indices
@deprecate safe_indices all_indices
@deprecate flatmatrix flatarray
@deprecate flatvector flatarray
@deprecate fastmatrix fastarray
@deprecate fastvector fastarray

@deprecate bcastdim(a::Integer, b::Integer) bcastsize(a, b)
@deprecate bcastdims bcastsize


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
