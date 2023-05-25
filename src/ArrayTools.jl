module ArrayTools

export
    ..,
    ArrayAxes,
    ArrayAxis,
    ArraySize,
    MaybeArrayAxes,
    MaybeArrayAxis,
    RubberIndex,
    UndefinedType,
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
    to_axes,
    to_axis,
    to_int,
    to_size,
    to_type,
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
    is_fast_array,
    # Macros:
    @assert_same_indices

using Base: OneTo, axes1, tail
import Base: dotview, getindex, setindex!, to_indices

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
