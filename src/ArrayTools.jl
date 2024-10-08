module ArrayTools

export
    ..,
    ArraySize,
    MaybeArrayAxes,
    MaybeArrayAxis,
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
    reversemap,
    all_indices,
    same_axes,
    same_size,
    same_standard_size,
    split_interval,
    standard_size,
    strictmap!,
    to_int,
    # to_type, # FIXME: deprecated, use `as(T,x)` instead

    # Re-exports from TypeUtils:
    ArrayAxes,
    ArrayAxis,
    ArrayShape,
    as,
    as_array_axes,
    as_array_axis,
    as_array_dim,
    as_array_size,
    as_eltype,
    convert_eltype,
    promote_eltype,

    # Storage trait:
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
    @assert_same_axes

using TypeUtils

using Base: OneTo, axes1, tail
import Base: dotview, getindex, setindex!, to_indices

@deprecate to_type(T, x) as(T, x) true

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
