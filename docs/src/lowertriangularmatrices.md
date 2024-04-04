# [LowerTriangularMatrices](@id lowertriangularmatrices)

LowerTriangularMatrices is a submodule that has been developed for SpeedyWeather.jl which is
technically independent (SpeedyWeather.jl however imports it and so does SpeedyTransforms)
and can also be used without running simulations. It is just not put into its own respective repository.

This module defines `LowerTriangularArray`, a lower triangular matrix format, which in contrast to
[`LinearAlgebra.LowerTriangular`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.LowerTriangular) does not store the entries above the diagonal. SpeedyWeather.jl
uses `LowerTriangularArray` which is defined as a subtype of `AbstractArray` to store
the spherical harmonic coefficients (see [Spectral packing](@ref)). For 2D `LowerTriangularArray` the alias `LowerTriangularMatrix` exists. Higher dimensional `LowerTriangularArray` are 'batches' of 2D `LowerTriangularMatrix`. So, for example a ``(10\times 10\times 10)`` `LowerTriangularArray` holds 10 `LowerTriangularMatrix` of size ``(10\times 10)`` in one array.  

## Creation of `LowerTriangularArray` 

A `LowerTriangularMatrix` and `LowerTriangularArray` can be created using `zeros`, `ones`, `rand`, or `randn`
```@repl LowerTriangularMatrices
using SpeedyWeather.LowerTriangularMatrices

L = rand(LowerTriangularMatrix{Float32}, 5, 5)
L2 = rand(LowerTriangularArray{Float32}, 5, 5, 5)
```
or the `undef` initializors `LowerTriangularMatrix{Float32}(undef, 3, 3)`, `.
The element type is arbitrary though, you can use any type `T` too.

Alternatively, it can be created through conversion from `Array`, which drops the upper triangle
entries and sets them to zero
```@repl LowerTriangularMatrices
M = rand(Float16, 3, 3)
L = LowerTriangularMatrix(M)

M2 = rand(Float16, 3, 3, 5)
L2 = LowerTriangularArray(M)
```

## Indexing `LowerTriangularArray`

`LowerTriangularArray` supports two types of indexing: 1) by denoting two indices, column and row `[l, m, ..]`
or 2) by denoting a single index `[lm, ..]`. The double index works as expected
```@repl LowerTriangularMatrices
L

L[2, 2]
```
But the single index skips the zero entries in the upper triangle, i.e.
```@repl LowerTriangularMatrices
L[4]
```
which, important, is different from single indices of an `AbstractMatrix`
```@repl LowerTriangularMatrices
Matrix(L)[4]
```
In performance-critical code a single index should be used, as this directly maps
to the index of the underlying data vector. The double index is somewhat slower
as it first has to be converted to the corresponding single index.

Consequently, many loops in SpeedyWeather.jl are build with the following structure
```@repl LowerTriangularMatrices
n, m = size(L)

ij = 0

for j in 1:m, i in j:n
    ij += 1
    L[ij] = i+j
end
```
which loops over all lower triangle entries of `L::LowerTriangularArray` and the single
index `ij` is simply counted up. However, one could also use `[i, j]` as indices in the
loop body or to perform any calculation (`i+j` here).
An iterator over all entries in the lower triangle can be created by
```@repl LowerTriangularMatrices
for ij in eachindex(L)
    # do something
end
```
The `setindex!` functionality of matrixes will throw a `BoundsError` when trying to write
into the upper triangle of a `LowerTriangularArray`, for example
```@repl LowerTriangularMatrices
L[2, 1] = 0    # valid index

L[1, 2] = 0    # invalid index in the upper triangle
```

Higher dimensional `LowerTriangularArray` can be indexed with multidimensional array indices 
like most other arrays types. Both the single index and the double index for the lower 
triangle work as shown here
```@repl LowerTriangularMatrices
L = rand(LowerTriangularMatrix{Float32}, 3, 3, 5)

L[2, 1] # second lower triangle element of the first lower triangle matrix 

L[2, 1, 1] # (2,1) element of the first lower triangle matrix 
```
The `setindex!` functionality follows accordingly. 

## Linear algebra with `LowerTriangularArray`

The [LowerTriangularMatrices](@ref lowertriangularmatrices) module's main purpose is not linear algebra, and it's
implementation may not be efficient, however, many operations work as expected
```@repl LowerTriangularMatrices
L = rand(LowerTriangularMatrix{Float32}, 3, 3)

L + L

L * L
```
Note, however, that the latter includes a conversion to `Matrix`, which is true for many
operations, including `inv` or `\`. Hence when trying to do more sophisticated linear
algebra with `LowerTriangularMatrix` we quickly leave lower triangular-land and go
back to normal matrix-land.

## GPU 

`LowerTriangularArray{T,N,ArrayType}` wraps around an array of type `ArrayType`. If this array is a GPU array (e.g. `CuArray`), all operations are performed on GPU as well. The implementation was written so that scalar indexing is avoided in almost all cases, so that GPU operation should be performant. To use `LowerTriangularArray` on GPU you can e.g. just `adapt` an exisiting `LowerTriangularArray`.
```julia 
L = rand(LowerTriangularArray{Float32}, 5, 5, 5)
L_gpu = adapt(CuArray, L)
```

## Function and type index

```@autodocs
Modules = [SpeedyWeather.LowerTriangularMatrices]
```
