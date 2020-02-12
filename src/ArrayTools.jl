module ArrayTools

export
    ..,
    RubberIndex,
    all_match,
    allof,
    anyof,
    bcastcopy,
    bcastdim,
    bcastdims,
    bcastlazy,
    cartesian_indices,
    checkdimensions,
    colons,
    Dimensions,
    dimensions,
    axis_limits,
    has_standard_indexing,
    noneof,
    promote_eltype,
    reversemap,
    safe_indices,
    same_axes,
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
@deprecate safeindices safe_indices
@deprecate flatmatrix flatarray
@deprecate flatvector flatarray
@deprecate fastmatrix fastarray
@deprecate fastvector fastarray

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
