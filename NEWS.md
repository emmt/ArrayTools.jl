# User visible changes in `ArrayTools` package

## Version 0.2.5

-  Relax signature of `reversemap`.

## Version 0.2.4

- Extend `promote_eltype` to any arguments extending the `eltype` method.
- Replace TravisCI by GitHub actions.

## Version 0.2.3

- New aliases `ArrayAxis` to `AbstractUnitRange{Int}` and `MaybeArrayAxis` to
  `AbstractUnitRange{Integer}` to represent an argument that is a valid array
  axis or eligible to be an array axis.  Similarly `ArrayAxes` and
  `MaybeArrayAxes` are aliases to tuples of `ArrayAxis` and `MaybeArrayAxis` to
  represent an argument that is a valid tuple of array axes or eligible to be a
  tuple of array axes.  New methods `to_axis` and `to_axes` are provided to
  respectively convert their argument(s) to instances of `ArrayAxis` and
  `ArrayAxes`.

- New method `to_type(T,x)` to convert argument `x` to type `T` with a type
  assertion on the result.

## Version 0.2.2

- Documentation has been improved.

- Some code ahs been simplified.


## Version 0.2.1

- `isfastarray` and `isflatarray` deprecated in favor of `is_fast_array` and
  `is_flat_array` which are more readable and thus less confusing.  `fastarray`
  and `flatarray` deprecated in favor of `to_fast_array` and `to_flat_array`
  which are more clear about their purpose.

- New method `same_size` to get the common size of arrays checking that they
  all have the same size.

- `same_dimensions` deprecated in favor of `same_standard_size` which is more
  explicit.

- `check_dimensions` deprecated in favor of `check_size` which yields the
  number of elements and throws an `ArgumentError` exception if not all
  dimensions are nonnegative.

- `dimensions(A::AbstractArray)` has been deprecated in favor of
  `standard_size(A)` while `dimensions` has been deprecated in favor of
  `to_size` for other types of arguments.  Union `Dimensions` has been
  deprecated favor of `ArraySize`.

- New method `to_int` for quick conversion of values to `Int`.

- `bcastdim` and `bcastdims` have been deprecated in favor of `bcastsize`.

- `checkdimensions` has been deprecated in favor of `check_dimensions`.

- Rubber index operator `…` has been renamed as `..` (which is easier to type) as in
  [`EllipsisNotation`](https://github.com/ChrisRackauckas/EllipsisNotation.jl).
