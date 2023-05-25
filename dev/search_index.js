var documenterSearchIndex = {"docs":
[{"location":"storage/#Array-storage","page":"Array storage","title":"Array storage","text":"","category":"section"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"When calling (with ccall) a compiled function coded in another language (C, FORTRAN, etc.), you have to make sure that array arguments have the same storage as assumed by these languages so that it is safe to pass the pointer of the array to the compiled function.","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"To check whether the elements of an array A are stored in memory contiguously and in column-major order, call:","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"is_flat_array(A) -> bool","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"which yield a bollean result.  Several arguments can be checked in a single call:","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"is_flat_array(A, B, C, ...)","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"is the same as:","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"is_flat_array(A) && is_fast_array(B) && is_fast_array(C) && ...","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"In order to get an array with such flat storage and possibly with a given element type T, call the to_flat_array method:","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"to_flat_array([T = eltype(A),] A)","category":"page"},{"location":"storage/","page":"Array storage","title":"Array storage","text":"which just returns A if the requirements hold for it, or which converts A to a suitable array form.","category":"page"},{"location":"broadcasting/#Broadcasting-of-arrays","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"","category":"section"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"ArrayTools provides the bcastlazy and bcastcopy methods for broadcasting arrays (that is expand an array to change the lengths of unit dimensions or to add trailing dimensions) with optional element type conversion.","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"For instance:","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"bcastlazy(A, [T=eltype(A),] dims...)","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"yields a flat array (see is_flat_array) of type T and dimensions dims whose values are given by A and according to type conversion and broadcasting rules (see broadcast standard method).","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"Making a copy of A is avoided by bcastlazy (hence its name), and A or a view on A is returned if A is already an array with the correct element type and dimensions or if it can be reshaped (by the reshape method) to match the constraints.  This means that the result may share the same contents as A.  If an array that does not share its contents with A is needed, then call:","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"bcastcopy(A, [T=eltype(A),] dims...)","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"instead of  bcastlazy.","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"The result returned by bcastlazy and  bcastcopy has 1-based indices and contiguous elements which is suitable for fast linear indexing.","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"The bcastsize method may be useful to determine dimensions when applying broadcasting rules:","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"bcastsize(size(A), size(B), ...) -> siz\nbcastsize(a, b) -> c","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"The first one yields the size siz of the array that would result from applying broadcasting rules to arguments A, B, etc.  The result is a tuple of integers (of type Int).  You may call check_size to ensure that the result is a valid list on nonnegative dimensions.","category":"page"},{"location":"broadcasting/","page":"Broadcasting of arrays","title":"Broadcasting of arrays","text":"The second applies broadcasting rules for a single dimension, throwing an exception if dimensions a and b are not compatible according to these rules.","category":"page"},{"location":"rubberindex/#Rubber-indices","page":"Rubber indices","title":"Rubber indices","text":"","category":"section"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"The constants .. and … (type \\dots and hit the tab key) can be used in array indexation to left or right justify the other indices.  For instance, assuming A is a 3×4×5×6 array, then all the following equalities hold:","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"A[..]          === A # the two are the same object\nA[..,3]         == A[:,:,:,3]\nA[2,..]         == A[2,:,:,:]\nA[..,2:4,5]     == A[:,:,2:4,5]\nA[:,2:3,..]     == A[:,2:3,:,:]\nA[2:3,..,1,2:4] == A[2:3,:,1,2:4]","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"As you can see, the advantage of the rubber index .. is that it automatically expands as the number of colons needed to have the correct number of indices.  The expressions are also more readable.  The idea comes from the Yorick language by Dave Munro.  Similar notation exists in NumPy.","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"The rubber index may also be used for setting values.  For instance:","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"A[..] .= 1         # to fill A with ones\nA[..,3] = A[..,2]  # to copy A[:,:,:,2] in A[:,:,:,3]\nA[..,3] .= A[..,2] # idem but faster\nA[2,..] = A[3,..]  # to copy A[3,:,:,:] in A[2,:,:,:]\nA[..,2:4,5] .= 7   # to set all elements in A[:,:,2:4,5] to 7","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"Leading/trailing indices may be specified as Cartesian indices (of type CartesianIndex).","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"Technically, the constant .. is defined as RubberIndex() where RubberIndex is the singleron type that represents any number of indices.","category":"page"},{"location":"rubberindex/","page":"Rubber indices","title":"Rubber indices","text":"Call colons(n) if you need a n-tuple of colons :.  When n is known at compile time, it is faster to call colons(Val(n)).","category":"page"},{"location":"install/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install/","page":"Installation","title":"Installation","text":"ArrayTools can be installed as any other offical Julia packages:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"… pkg> add ArrayTools","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"where … pkg> represents the package manager prompt (the ellipsis … denote your current environment).  To start Julia's package manager, launch Julia and, at the REPL of Julia, hit the ] key; you should get the above … pkg> prompt.  To revert to Julia's REPL, just hit the Backspace key at the … pkg> prompt.","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"To check whether the ArrayTools package works correctly, type:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"… pkg> test ArrayTools","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"Later, to update to the last version (and run tests), you can type:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"… pkg> update ArrayTools\n… pkg> test ArrayTools","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"If something goes wrong, it may be because you already have an old version of ArrayTools.  Uninstall ArrayTools as follows:","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"… pkg> rm ArrayTools\n… pkg> gc\n… pkg> add https://github.com/emmt/ArrayTools.jl","category":"page"},{"location":"install/","page":"Installation","title":"Installation","text":"before re-installing.","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"The following provides detailled documentation about types and methods provided by the ArrayTools package.  This information is also available from the REPL by typing ? followed by the name of a method or a type.","category":"page"},{"location":"reference/#Broadcasting","page":"Reference","title":"Broadcasting","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"bcastlazy\nbcastcopy\nbcastsize","category":"page"},{"location":"reference/#ArrayTools.bcastlazy","page":"Reference","title":"ArrayTools.bcastlazy","text":"bcastlazy(A, [T=eltype(A),] dims...)\n\nyields a flat array of type T and dimensions dims whose values are given by A according to type conversion and broadcasting rules (see broadcast method). Compared to bcastcopy, making a copy of A is avoided if it is already an array with the correct type of elements and dimensions or if it can be reshaped (by the reshape method) to the correct type and dimensions. This means that the result may share the same contents as A. Argument A can be a scalar or an array with 1-based indices. The result has 1-based indices and contiguous elements which is suitable for fast linear indexing.\n\nSee also bcastcopy, bcastsize.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.bcastcopy","page":"Reference","title":"ArrayTools.bcastcopy","text":"bcastcopy(A, [T=eltype(A),] dims...)\n\nyields a new array of element type T and dimensions dims whose values are given by A according to type conversion and broadcasting rules (like for the broadcast method). Compared to bcastlazy, it is guaranteed that the returned array does not share its contents with A.\n\nArgument A can be a scalar value or an array.\n\nSee also bcastlazy, bcastsize.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.bcastsize","page":"Reference","title":"ArrayTools.bcastsize","text":"bcastsize(size(A), size(B), ...) -> siz\n\nyields the size siz of the array that would result from applying broadcasting rules (see broadcast method) to arguments A, B, etc. The result is a tuple of integers (of type Int). Call check_size if you want to also make sure that the result is a list of valid dimensions.\n\nThe method can also be applied to a single dimension:\n\nbcastsize(a, b) -> c\n\nto yield the dimension c gievn by broadcasting dimensions a and b throwing an exception if dimensions are not compatible according to broadcasting rules. This is the same as Base.Broadcasting._bcs1 but it takes care of converting to Int.\n\nSee also standard_size, check_size, bcastcopy, bcastlazy.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Indexing","page":"Reference","title":"Indexing","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"ArraySize\nRubberIndex\ncolons\nIndexingType\nFastIndexing\nAnyIndexing\nis_fast_array\nto_fast_array\nto_size\ncheck_size\nsame_size\nsame_standard_size\nstandard_size\nsame_axes\n@assert_same_indices\nall_indices\ncartesian_indices\ncommon_indices\nhas_standard_indexing","category":"page"},{"location":"reference/#ArrayTools.ArraySize","page":"Reference","title":"ArrayTools.ArraySize","text":"ArraySize\n\nis the union of types eligible to define array size. Calling [to_size](@ref)(dims) on any argument dims such that isa(dims,ArraySize) is true yields an array size in canonical form, that is an instance of Dims{N} which is an alias for an N-tuple of Int.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.RubberIndex","page":"Reference","title":"ArrayTools.RubberIndex","text":"RubberIndex is the singleron type that represents any number of indices. The constant .. is defined as RubberIndex() and can be used in array indexation to left and/or right justify the other indices. For instance, assuming A is a 3×4×5×6 array, then all the following equalities hold:\n\nA[..]           == A[:,:,:,:]\nA[..,3]         == A[:,:,:,3]\nA[2,..]         == A[2,:,:,:]\nA[..,2:4,5]     == A[:,:,2:4,5]\nA[2:3,..,1,2:4] == A[2:3,:,1,2:4]\n\nAs you can see, the advantage of the rubber index .. is that it automatically expands as the number of colons needed to have the correct number of indices. The expressions are also more readable.\n\nThe rubber index may also be used for setting values. For instance:\n\nA[..] .= 1         # to fill A with ones\nA[..,3] = A[..,2]  # to copy A[:,:,:,2] in A[:,:,:,3]\nA[..,3] .= A[..,2] # idem but faster\nA[2,..] = A[3,..]  # to copy A[3,:,:,:] in A[2,:,:,:]\nA[..,2:4,5] .= 7   # to set all elements in A[:,:,2:4,5] to 7\n\nLeading/trailing indices may be specified as Cartesian indices (of type CartesianIndex).\n\nwarning: Warning\nThere are two known limitations:The end reserved word can only be used in intervals specified before the rubber index but not after. This limitation is due to the Julia parser cannot be avoided.\nAt most 9 indices can be specified before the rubber index. This can be extended by editing the source code.\n\nSee also: colons.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.colons","page":"Reference","title":"ArrayTools.colons","text":"colons(n)\n\nyields a n-tuple of colons : (a.k.a. Colon()).\n\nWhen n is known at compile time, it is faster to call:\n\ncolons(Val(n))\n\nThis method is suitable to extract sub-arrays of build views when some kind of rubber index is needed. For instance:\n\nslice(A::AbstractArray{T,N}, i::Integer) where {T,N} =\n    A[colons(Val(N-1))..., i]\n\ndefines a function that returns the i-th slice of A assuming index i refers the last index of A. Using the rubber-index .., a shorter definition is:\n\nslice(A::AbstractArray, i) = A[.., i]\n\nwhich is also able to deal with multiple trailing indices if i is a CartesianIndex.\n\nSee also: .., RubberIndex.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.IndexingType","page":"Reference","title":"ArrayTools.IndexingType","text":"IndexingType(A)\n\nyields one of the singletons FastIndexing() or AnyIndexing() to indicate whether or not array A has standard 1-based indices and can be efficiently indexed by one integer (even if A is multidimensional) and column-major ordering is used to access the elements of A.\n\nThis method can be extended for custom array types to quickly return the correct answer.\n\nSee also is_fast_array, to_fast_array.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.FastIndexing","page":"Reference","title":"ArrayTools.FastIndexing","text":"FastIndexing()\n\nyields the indexing type of fast arrays. See IndexingType.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.AnyIndexing","page":"Reference","title":"ArrayTools.AnyIndexing","text":"AnyIndexing()\n\nyields the indexing type of non-fast arrays. See IndexingType.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.is_fast_array","page":"Reference","title":"ArrayTools.is_fast_array","text":"is_fast_array(A)\n\nyields whether array A has standard 1-based indices and is efficiently indexed by linear indices.\n\nSeveral arguments can be checked in a single call:\n\nis_fast_array(A, B, C, ...)\n\nis the same as:\n\nis_fast_array(A) && is_fast_array(B) && is_fast_array(C) && ...\n\nSee also IndexingType, to_fast_array, is_flat_array.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.to_fast_array","page":"Reference","title":"ArrayTools.to_fast_array","text":"to_fast_array([T=eltype(A),] A)\n\nlazily yields a fast array equivalent to A with element type T. A fast array has standard 1-based indices and is efficiently indexed by linear indices. If A is already a fast array with element type T, A is returned; otherwise, A is converted into an Array which is returned.\n\nSee also is_fast_array, IndexingType, to_flat_array.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.to_size","page":"Reference","title":"ArrayTools.to_size","text":"to_size(dims)\n\nconverts dims to an instance of Dims{N} which is an alias for an N-tuple of Int. Argument dims can be a scalar integer or a tuple of integers. Argument dims is returned if already of the correct type. This method may also be called as:\n\nto_size(dim1, dim2, ...)\n\nto let to_size deal with a variable number of arguments.\n\nThe union ArraySize matches the types of acceptable argument(s) for to_size(arg): scalar integers and tuples of integers.\n\nThis method is intended for fast conversion, call check_size(dims) to verify that all dimensions in dims are nonnegative.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.check_size","page":"Reference","title":"ArrayTools.check_size","text":" check_size(siz) -> len\n\nchecks the validity of the array size siz and yields the corresponding number of elements (throwing an ArgumentError exception if this is not the case). To be a valid array size, the values of siz must all be nonnegative.\n\nSee also to_size.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.same_size","page":"Reference","title":"ArrayTools.same_size","text":"same_size(A, B...) -> size(A)\n\nchecks whether arrays A, B, etc., all have the same size which is returned. A DimensionMismatch exception is thrown if array sizes are not all identical.\n\nSee also same_standard_size, same_axes.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.same_standard_size","page":"Reference","title":"ArrayTools.same_standard_size","text":"same_standard_size(A, B...) -> size(A)\n\nchecks whether arrays A, B, etc., all have standard indexing and the same size which is returned. If array sizes are not all identical, a DimensionMismatch exception is thrown. If arrays have non-standard indexing (that is indices not starting at index one), an ArgumentError exception is thrown.\n\nSee also standard_size, has_standard_indexing.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.standard_size","page":"Reference","title":"ArrayTools.standard_size","text":"standard_size(A) -> size(A)\n\nyields the list of dimensions of A, that is size(A), throwing an ArgmentError exception if A does not have standard 1-based indices.\n\nSee also has_standard_indexing, same_standard_size.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.same_axes","page":"Reference","title":"ArrayTools.same_axes","text":"same_axes(A, B...) -> axes(A)\n\nchecks whether arrays A, B, etc., have the same axes and returns them. A DimensionMismatch exception is thrown if axes are not all identical.\n\nSee also same_size, all_indices.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.@assert_same_indices","page":"Reference","title":"ArrayTools.@assert_same_indices","text":"@assert_same_indices A B ...\n\nthrows a DimensionMismatch exception if arrays A, B, etc. do not have the same indices.\n\n\n\n\n\n","category":"macro"},{"location":"reference/#ArrayTools.all_indices","page":"Reference","title":"ArrayTools.all_indices","text":"all_indices(A...)\n\nyields an iterable object for visiting each index of array(s) A in an efficient manner. For array types that have opted into fast linear indexing (like Array), this is simply the range 1:length(A). For other array types, return a specialized Cartesian range to efficiently index into the array(s) with indices specified for every dimension.\n\nIf more than one AbstractArray argument are supplied, all_indices will create an iterable object that is fast for all arguments (a UnitRange if all inputs have fast linear indexing, a CartesianIndices otherwise). A DimensionMismatch exception is thrown if the arrays have different axes so that it is always safe to use @inbounds in front of a loop like:\n\nfor i in all_indices(A, B, C, D)        A[i] = B[i]*C[i] + D[i]    end\n\nwhen A, B etc. are all (abstract) arrays.\n\nThis method is similar to eachindex except that a DimensionMismatch exception is thrown if arrays have different axes. For linearly indexed arrays, eachindex only checks that they have the same linear index range (that is the same number of elements, not the same shape).\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.cartesian_indices","page":"Reference","title":"ArrayTools.cartesian_indices","text":"cartesian_indices(A)\ncartesian_indices((n1, n2, ...))\ncartesian_indices((i1:j1, i2:j2, ...))\ncartesian_indices(CartesianIndex(i1, i2, ...), CartesianIndex(j1, j2, ...))\ncartesian_indices(R)\n\nall yield an instance of CartesianIndices suitable for multi-dimensional indexing of respectively: all the indices of array A, a multi-dimensional array of dimensions (n1,n2,...), a multi-dimensional region whose first and last indices are (i1,i2,...) and (j1,j2,...) or a Cartesian region defined by R, an instance of CartesianIndices.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.common_indices","page":"Reference","title":"ArrayTools.common_indices","text":"Assuming A and B are arrays with N dimensions:\n\ncommon_indices(A, B) -> inds\n\nyields the set of all the indices that are valid for both A and B. The result is similar to axes(A) or axes(B), that is an N-tuple of integer valued unit ranges.\n\nAn offset k with a sign may be specified:\n\ncommon_indices(A, B, ±, k)\n\nto obtain the set of all indices i such that A[i] and B[i ± k] are valid and where here and above ± is either + or -. Offset k can be a tuple of integers or a Cartesian index.\n\nArguments A and B may be both tuples of indices or index ranges or both scalar or index range which specify the size or the axes of the arrays to be indexed. This is used in the following example, where we want to do A[i] = B[i]*C[i + k] given the offset k and for all valid indices i:\n\nI = common_indices(same_axes(A, B), axes(C), +, k)\n@inbounds @simd for i in CartesianIndices(I)\n   A[i] = B[i]*C[i + k]\nend\n\nNote that same_axes(A,B) is called to get the axes of A and B while asserting that they are the same, as a result no bound checking is necessary and the loop can be optimized for vectorization.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.has_standard_indexing","page":"Reference","title":"ArrayTools.has_standard_indexing","text":"has_standard_indexing(A)\n\nreturn true if the indices of A start with 1 along all axes. Can be called with multiple arguments:\n\nhas_standard_indexing(A, B, ...)\n\nis equivalent to:\n\nhas_standard_indexing(A) && has_standard_indexing(B) && ...\n\nOpposite of Base.has_offset_axes which is not available in version of Julia older than 0.7.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Storage","page":"Reference","title":"Storage","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"StorageType\nAnyStorage\nFlatStorage\nis_flat_array\nto_flat_array","category":"page"},{"location":"reference/#ArrayTools.StorageType","page":"Reference","title":"ArrayTools.StorageType","text":"StorageType(A)\n\nyields the type of storage of the elements of argument A. If A is a flat array, that is an array with contiguous elements in column-major order and first element at index 1, the singleton FlatStorage() is returned; otherwise, the singleton AnyStorage() is returned.\n\nThis method can be extended for custom array types to quickly return the correct answer.\n\nSee also is_flat_array, to_flat_array.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.AnyStorage","page":"Reference","title":"ArrayTools.AnyStorage","text":"AnyStorage()\n\nyields the storage type of a non-flat arrays. See StorageType.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.FlatStorage","page":"Reference","title":"ArrayTools.FlatStorage","text":"FlatStorage()\n\nyields the storage type of flat arrays. See StorageType.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ArrayTools.is_flat_array","page":"Reference","title":"ArrayTools.is_flat_array","text":"is_flat_array(A) -> boolean\n\nyields whether array A can be indexed as a flat array, that is an array with contiguous elements in column-major order and first element at index 1. This also means that A has 1-based indices along all its dimensions.\n\nSeveral arguments can be checked in a single call:\n\nis_flat_array(A, B, C, ...)\n\nis the same as:\n\nis_flat_array(A) && is_flat_array(B) && is_flat_array(C) && ...\n\nSee also StorageType, to_flat_array, is_fast_array, has_standard_indexing.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.to_flat_array","page":"Reference","title":"ArrayTools.to_flat_array","text":"to_flat_array([T=eltype(A),] A)\n\nlazily yields a flat array based on A, that is an array with contiguous elements in column-major order and first element at index 1. Optional argument T is to specify the element type of the result. Argument A is returned if it is already a flat array with the requested element type; otherwise, convert method is called to produce the result (an Array{T} in that case).\n\nSee also is_flat_array, to_fast_array.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Pseudo-arrays","page":"Reference","title":"Pseudo-arrays","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"ArrayTools.PseudoArray","category":"page"},{"location":"reference/#ArrayTools.PseudoArrays.PseudoArray","page":"Reference","title":"ArrayTools.PseudoArrays.PseudoArray","text":"Abstract type PseudoArray{T,N,S} is to be derived by types that want to provide an array-like interface. Parameter T is the element type, parameter N is the number of dimensions and parameter S is the index style: IndexCartesian or IndexLinear.\n\nnote: Note\nThe indexing style must be part of the signature because it must be possible to call IndexStyle() on the data type not the instance. Another possibility would have been to have the type of the embedded array be part of the signature but this is more restrictive.\n\nAlias LinearArray{T,N} is an abstract type that can be derived by types that want to provide an array-like interface with array values stored in an array whose index style is linear.\n\nUsage can be as simple as:\n\nstruct CustomArray{T,N,...} <: LinearArray{T,N}\n    arr::Array{T,N} # can be any array type with linear index style\n    ...             # anything else\nend\n\n@inline Base.parent(A::CustomArray) = A.arr\n\nAs a result, instances of CustomArray{T,N} will be seen as instances of AbstractArray{T,N} and behave as if they implement linear indexing. Apart from the needs to extend the Base.parent method, the interface to LinearArray{T,N} should provide any necessary methods for indexation, getting the dimensions, the element type, etc. for the derived custom type. You may however override these definitions by more optimized or more suitable methods specialized for your custom array-like type.\n\nSimilarly, alias CartesianArray{T,N} is an abstract type that can be derived by types that want to provide an array-like interface with array values stored in an array whose index style is Cartesian. For such array-like object, index checking requires an efficient implementation of the Base.axes() method which you may have to specialize. The default implementation is:\n\n@inline Base.axes(A::PseudoArray) = axes(parent(A))\n\n\n\n\n\n","category":"type"},{"location":"reference/#Utilities","page":"Reference","title":"Utilities","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"all_match\nallof\nanyof\naxis_limits\nnoneof\npromote_eltype\nreversemap\nsplit_interval\nstrictmap!","category":"page"},{"location":"reference/#ArrayTools.all_match","page":"Reference","title":"ArrayTools.all_match","text":"all_match(val, f, args...) -> bool\n\nyields as soon as possible (short-circuit) whether f(arg) == val for each argument arg in args.... The returned value is true if there are no arguments after f.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.allof","page":"Reference","title":"ArrayTools.allof","text":"allof(f, args...) -> Bool\n\nchecks whether predicate function f returns true for all arguments in args..., returning false as soon as possible (short-circuiting).\n\nallof(args...) -> Bool\n\nchecks whether all arguments args... are true, returning false as soon as possible (short-circuiting). Arguments can be booleans or arrays of booleans. The latter are considered as true if all their elements are true and are considered as false otherwise (if any of their elements are false). Arguments can also be iterables to check whether all their values are true. An empty iterable is considered as true.\n\nThis method can be much faster than all(f, args) or all(args) because its result may be determined at compile time. However, missing values are not considered as special.\n\nSee also all, anyof, noneof.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.anyof","page":"Reference","title":"ArrayTools.anyof","text":"anyof(f, args...) -> bool\n\nchecks whether predicate function f returns true for any argument args..., returning true as soon as possible (short-circuiting).\n\nanyof(args...) -> bool\n\nchecks whether all arguments args... are true, returning false as soon as possible (short-circuiting). Arguments can be booleans or arrays of booleans. The latter are considered as true if any of their elements are true and are considered as false otherwise (if all their elements are false). Arguments can also be iterables to check whether any of their values are true. An empty iterable is considered as false.\n\nThis method can be much faster than any(f, args) or any(args) because its result may be determined at compile time. However, missing values are not considered as special.\n\nSee also any, allof, noneof.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.axis_limits","page":"Reference","title":"ArrayTools.axis_limits","text":"axis_limits(I) = (i0,i1)\n\nyields the limits i0 and i1 of index range I as a 2-tuple of Int's and such that i0:i1 represents the same indices as I (although not in the same order if step(I) < 0). If step(I) is not equal to ±1, an ArgumentError exception is thrown.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.noneof","page":"Reference","title":"ArrayTools.noneof","text":"noneof(f, args...) -> bool\n\nchecks whether predicate f returns false for all argument args..., while\n\nnoneof(args...) -> bool\n\nchecks whether all argument args... are false.\n\nSee also any, allof, noneof.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.promote_eltype","page":"Reference","title":"ArrayTools.promote_eltype","text":"promote_eltype(args...)\n\nyields the promoted element type of its arguments. Arguments args... may be anything implementing the eltype method.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.reversemap","page":"Reference","title":"ArrayTools.reversemap","text":"reversemap(f, args)\n\napplies the function f to arguments args in reverse order and return the result. For now, the arguments args must be in the form of a simple tuple and the result is the tuple: (f(args[end]),f(args[end-1]),...,f(args[1]).\n\nAlso see: map, ntuple.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.split_interval","page":"Reference","title":"ArrayTools.split_interval","text":"split_interval(I, J, +/-, k) -> Ia, Ib, Ic\n\ngiven unit ranges I and J and offset ±k, yields 3 unit ranges, such that Ia ∪ Ib ∪ Ic = I and:\n\n∀ i ∈ Ia, i ± k < first(J);\n∀ i ∈ Ib, i ± k ∈ J;\n∀ i ∈ Ic, i ± k > last(J).\n\nUnit ranges may be replaced by their first and last values:\n\nsplit_interval(first(I), last(I), first(J), last(J), +/-, k)\n\nyields the same result as above.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ArrayTools.strictmap!","page":"Reference","title":"ArrayTools.strictmap!","text":"strictmap!(dst, f, src) -> dst\n\ndoes dst[i] = f(src[i]) for all indices i and returns dst. Arguments dst and src must have the same axes.\n\nExcept for the strict condition on the axes, this method is similar to map!(f,dst,src).\n\n\n\n\n\n","category":"function"},{"location":"indexing/#Array-Indexing","page":"Array Indexing","title":"Array Indexing","text":"","category":"section"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"The all_indices method takes any number of array arguments and yields an efficient iterator for visiting all indices each index of the arguments. Its behavior is similar to that of eachindex method except that all_indices throws a DimensionMismatch exception if linearly indexed arrays have different axes while eachindex just checks that such arrays have the same number of elements. It is always safe to specify @inbounds (and @simd or @turbo from the LoopVectorization package) for a loop like:","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"for i in all_indices(A, B, C, D)\n   A[i] = B[i]*C[i] + D[i]\nend","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"An alternative is to call the @assert_same_indices macro which throws a DimensionMismatch exception if the provided arguments are not arrays with the same indices. For example:","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"@assert_same_indices A B C D\n@inbounds for i in eachindex(A, B, C, D)\n   A[i] = B[i]*C[i] + D[i]\nend","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"The eachindex and all_indices methods are very useful when writing loops over array elements so as to be agnostic to which specific indexing rule is the most suitable. Some algorithms are however more efficient or easier to write if all involved arrays are indexed by a single 1-based index. ArrayTools provides is_fast_array to check whether arrays are suitable for fast indexing:","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"is_fast_array(A) -> bool","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"to check array A or:","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"is_fast_array(A, B, ...) -> bool","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"to check several arrays A, B, ... at the same time. ArrayTools also provides to_fast_array to convert an array to fast linear indexing if this is needed:","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"to_fast_array(A)","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"checks whether array A is suitable for fast indexing (by a single integer starting at 1); if it does, A is returned to the caller; otherwise, the contents of A is converted to a suitable array type implementing fast indexing and is returned to the caller.","category":"page"},{"location":"indexing/","page":"Array Indexing","title":"Array Indexing","text":"The method has_standard_indexing can be called to check whether one or several arrays have standard 1-based indexing though they not necessarily implement fast linear indexing.","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"ArrayTools is a Julia package which provides a number of methods and types to deal with the variety of array types (sub-types of AbstractArray) that exist in Julia and to help building custom array-like types without sacrificing performances.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"ArrayTools is hosted on GitHub.","category":"page"},{"location":"#Table-of-contents","page":"Introduction","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\"install.md\", \"storage.md\", \"indexing.md\", \"rubberindex.md\",\n         \"broadcasting.md\", \"reference.md\"]","category":"page"},{"location":"#Index","page":"Introduction","title":"Index","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"}]
}
