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
