# Rubber indices

The constants `..` and `…` (type `\dots` and hit the `tab` key) can be used in
array indexation to left or right justify the other indices. For instance,
assuming `A` is a `3×4×5×6` array, then all the following equalities hold:

```julia
A[..]          === A # the two are the same object
A[..,3]         == A[:,:,:,3]
A[2,..]         == A[2,:,:,:]
A[..,2:4,5]     == A[:,:,2:4,5]
A[:,2:3,..]     == A[:,2:3,:,:]
A[2:3,..,1,2:4] == A[2:3,:,1,2:4]
```

As you can see, the advantage of the *rubber index* `..` is that it
automatically expands as the number of colons needed to have the correct number
of indices. The expressions are also more readable. The idea comes from the
[`Yorick`](http://github.com/LLNL/yorick/) language by Dave Munro. Similar
notation exists in
[`NumPy`](https://numpy.org/doc/stable/user/basics.indexing.html).

The rubber index may also be used for setting values. For instance:

```julia
A[..] .= 1         # to fill A with ones
A[..,3] = A[..,2]  # to copy A[:,:,:,2] in A[:,:,:,3]
A[..,3] .= A[..,2] # idem but faster
A[2,..] = A[3,..]  # to copy A[3,:,:,:] in A[2,:,:,:]
A[..,2:4,5] .= 7   # to set all elements in A[:,:,2:4,5] to 7
```

Leading/trailing indices may be specified as Cartesian indices (of type
`CartesianIndex`).

Technically, the constant `..` is defined as `RubberIndex()` where
`RubberIndex` is the singleron type that represents any number of indices.

Call `colons(n)` if you need a `n`-tuple of colons `:`. When `n` is known at
compile time, it is faster to call `colons(Val(n))`.
