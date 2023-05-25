# Reference

The following provides detailed documentation about types and methods provided
by the `ArrayTools` package. This information is also available from the REPL
by typing `?` followed by the name of a method or a type.


## Broadcasting

```@docs
bcastlazy
bcastcopy
bcastsize
```

## Indexing

```@docs
ArraySize
RubberIndex
colons
IndexingType
FastIndexing
AnyIndexing
is_fast_array
to_fast_array
to_size
check_size
same_size
same_standard_size
standard_size
same_axes
@assert_same_axes
all_indices
cartesian_indices
common_indices
has_standard_indexing
```

## Storage

```@docs
StorageType
AnyStorage
FlatStorage
is_flat_array
to_flat_array
```

## Pseudo-arrays

```@docs
ArrayTools.PseudoArray
```


## Utilities

```@docs
all_match
allof
anyof
axis_limits
noneof
promote_eltype
reversemap
split_interval
strictmap!
```
