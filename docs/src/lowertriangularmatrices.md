# [LowerTriangularMatrices](@id lowertriangularmatrices)

LowerTriangularMatrices is a submodule that has been developed for SpeedyWeather.jl which is
technically independent (SpeedyWeather.jl however imports it and so does SpeedyTransforms)
and can also be used without running simulations. It is just not put into its own respective repository.

This module defines `LowerTriangularMatrix`, a lower triangular matrix, which in contrast to
[`LinearAlgebra.LowerTriangular`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.LowerTriangular) does not store the entries above the diagonal. SpeedyWeather.jl
uses `LowerTriangularMatrix` which is defined as a subtype of `AbstractMatrix` to store
the spherical harmonic coefficients (see [Spectral packing](@ref)). 

## Creation of `LowerTriangularMatrix`

A `LowerTriangularMatrix` can be created using `zeros`,`ones`,`rand`, or `randn`
```julia
julia> using SpeedyWeather.LowerTriangularMatrices

julia> L = rand(LowerTriangularMatrix{Float32},5,5)
5×5 LowerTriangularMatrix{Float32}:
 0.912744   0.0        0.0       0.0       0.0
 0.0737592  0.230592   0.0       0.0       0.0
 0.799679   0.0765255  0.888098  0.0       0.0
 0.670835   0.997938   0.505276  0.492966  0.0
 0.949321   0.193692   0.793623  0.152817  0.357968
```
or the `undef` initializor `LowerTriangularMatrix{Float32}(undef,3,3)`.
The element type is arbitrary though, you can use any type `T` too.

Alternatively, it can be created through conversion from `Matrix`, which drops the upper triangle
entries and sets them to zero
```julia
julia> M = rand(Float16,3,3)
3×3 Matrix{Float16}:
 0.2222  0.694    0.3452
 0.2158  0.04443  0.274
 0.9746  0.793    0.6294

julia> LowerTriangularMatrix(M)
3×3 LowerTriangularMatrix{Float16}:
 0.2222  0.0      0.0
 0.2158  0.04443  0.0
 0.9746  0.793    0.6294
```

## Indexing `LowerTriangularMatrix`

`LowerTriangularMatrix` supports two types of indexing: 1) by denoting two indices, column and row `[l,m]`
or 2) by denoting a single index `[lm]`. The double index works as expected
```julia
julia> L
3×3 LowerTriangularMatrix{Float16}:
 0.1499  0.0    0.0
 0.1177  0.478  0.0
 0.1709  0.756  0.3223

julia> L[2,2]
Float16(0.478)
```
But the single index skips the zero entries in the upper triangle, i.e.
```julia
julia> L[4]
Float16(0.478)
```
which, important, is different from single indices of an `AbstractMatrix`
```julia
julia> Matrix(L)[4]
Float16(0.0)
```
In performance-critical code a single index should be used, as this directly maps
to the index of the underlying data vector. The double index is somewhat slower
as it first has to be converted to the corresponding single index.

Consequently, many loops in SpeedyWeather.jl are build with the following structure
```julia
n,m = size(L)
ij = 0
for j in 1:m
    for i in j:n
        ij += 1
        L[ij] = i+j
    end
end
```
which loops over all lower triangle entries of `L::LowerTriangularMatrix` and the single
index `ij` is simply counted up. However, one could also use `[i,j]` as indices in the
loop body or to perform any calculation (`i+j` here).
An iterator over all entries in the lower triangle can be created by
```julia
for ij in eachindex(L)
    # do something
end
```
The `setindex!` functionality of matrixes will throw a `BoundsError` when trying to write
into the upper triangle of a `LowerTriangularMatrix`, for example
```julia
julia> L[2,1] = 0    # valid index
0

julia> L[1,2] = 0    # invalid index in the upper triangle
ERROR: BoundsError: attempt to access 3×3 LowerTriangularMatrix{Float32} at index [1, 2]
```

## Linear algebra with `LowerTriangularMatrix`

The [LowerTriangularMatrices](@ref lowertriangularmatrices) module's main purpose is not linear algebra, and it's
implementation may not be efficient, however, many operations work as expected
```julia
julia> L = rand(LowerTriangularMatrix{Float32},3,3)
3×3 LowerTriangularMatrix{Float32}:
 0.57649   0.0       0.0
 0.348685  0.875371  0.0
 0.881923  0.850552  0.998306

julia> L + L
3×3 LowerTriangularMatrix{Float32}:
 1.15298   0.0      0.0
 0.697371  1.75074  0.0
 1.76385   1.7011   1.99661

julia> L * L
3×3 Matrix{Float32}:
 0.332341  0.0       0.0
 0.506243  0.766275  0.0
 1.68542   1.59366   0.996616
```
Note, however, that the latter includes a conversion to `Matrix`, which is true for many
operations, including `inv` or `\`. Hence when trying to do more sophisticated linear
algebra with `LowerTriangularMatrix` we quickly leave lower triangular-land and go
back to normal matrix-land.

## Function and type index

```@autodocs
Modules = [SpeedyWeather.LowerTriangularMatrices]
```