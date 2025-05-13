# [LowerTriangularArrays](@id lowertriangularmatrices)

LowerTriangularArrays is a submodule that has been developed for SpeedyWeather.jl which is
technically independent (SpeedyWeather.jl however imports it and so does SpeedyTransforms)
and can also be used without running simulations. It is just not put into its own respective repository.

This module defines `LowerTriangularArray`, a lower triangular matrix format, which in contrast to
[`LinearAlgebra.LowerTriangular`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.LowerTriangular)
does not store the entries above the diagonal.
SpeedyWeather.jl uses `LowerTriangularArray` which is defined as a subtype of `AbstractArray` to store
the spherical harmonic coefficients (see [Spectral packing](@ref)).
For 2D `LowerTriangularArray` the alias `LowerTriangularMatrix` exists.
Higher dimensional `LowerTriangularArray` are 'batches' of 2D `LowerTriangularMatrix`.
So, for example a ``(10\times 10\times 10)`` `LowerTriangularArray` holds 10 `LowerTriangularMatrix` of size ``(10\times 10)`` in one array.  

!!! warn "LowerTriangularMatrix is actually a vector"
    `LowerTriangularMatrix` and `LowerTriangularArray` can in many ways be used very much like a `Matrix` or `Array`, however,
    because they unravel the lower triangle into a vector their dimensionality is one less than their `Array` counterparts.
    A `LowerTriangularMatrix` should therefore be treated as a vector rather than a matrix with some (limited) added
    functionality to allow for matrix-indexing (vector or flat-indexing is the default though). More details below.

## Creation of `LowerTriangularArray` 

A `LowerTriangularMatrix` and `LowerTriangularArray` can be created using `zeros`, `ones`, `rand`, or `randn`
```@repl LowerTriangularArrays
using SpeedyWeather.LowerTriangularArrays

L = rand(LowerTriangularMatrix{Float32}, 5, 5)
L2 = rand(LowerTriangularArray{Float32}, 5, 5, 5)
```
or the `undef` initializer `LowerTriangularMatrix{Float32}(undef, 3, 3)`.
The element type is arbitrary though, you can use any type `T` too.

Note how for a matrix both the upper triangle and the lower triangle are shown in the terminal.
The zeros are evident. However, for higher dimensional `LowerTriangularArray` we fall back to
show the unravelled first two dimensions. Hence, here, the first column is the first
matrix with 15 elements forming a 5x5 matrix, but the zeros are not shown.

Alternatively, it can be created through conversion from `Array`, which drops the upper triangle
entries and sets them to zero (which are not stored however)
```@repl LowerTriangularArrays
M = rand(Float16, 3, 3)
L = LowerTriangularMatrix(M)

M2 = rand(Float16, 3, 3, 2)
L2 = LowerTriangularArray(M2)
```

## Size of `LowerTriangularArray`

There are three different ways to describe the size of a `LowerTriangularArray`. For example with `L`

```@repl LowerTriangularArrays
L = rand(LowerTriangularMatrix, 5, 5)
```

we have (additional dimensions follow naturally thereafter)

### 1-based vector indexing (default)

```@repl LowerTriangularArrays
size(L)    # equivalently size(L, OneBased, as=Vector)
```

The lower triangle is unravelled hence the number of elements in the lower triangle is returned.

### 1-based matrix indexing

```@repl LowerTriangularArrays
size(L, as=Matrix)    # equivalently size(L, OneBased, as=Matrix)
```

If you think of a `LowerTriangularMatrix` as a matrix this is the most intuitive size of `L`, which,
however, does not agree with the size of the underlying data array (hence it is not the default).

### 0-based matrix indexing

Because `LowerTriangularArray`s are used to represent the coefficients of spherical harmonics which
are commonly indexed based on zero (i.e. starting with the zero mode representing the mean),
we also add `ZeroBased` to get the corresponding size.

```@repl LowerTriangularArrays
size(L, ZeroBased, as=Matrix)
```
which is convenient if you want to know the maximum degree and order of the spherical harmonics in `L`.
0-based vector indexing is not implemented.


## Indexing `LowerTriangularArray`

We illustrate the two types of indexing `LowerTriangularArray` supports.

- Matrix indexing, by denoting two indices, column and row `[l, m, ..]`
- Vector/flat indexing, by denoting a single index `[lm, ..]`.

The matrix index works as expected

```@repl LowerTriangularArrays
L

L[2, 2]
```
But the single index skips the zero entries in the upper triangle, i.e. a `2, 2` index points to the
same element as the index `6`
```@repl LowerTriangularArrays
L[6]
```
which, important, is different from single indices of an `AbstractMatrix`
```@repl LowerTriangularArrays
Matrix(L)[6]
```
which would point to the first element in the upper triangle (hence zero).

In performance-critical code a single index should be used, as this directly maps
to the index of the underlying data vector. The matrix index is somewhat slower
as it first has to be converted to the corresponding single index. 

Consequently, many loops in SpeedyWeather.jl are build with the following structure
```@repl LowerTriangularArrays
n, m = size(L, as=Matrix)

ij = 0
for j in 1:m, i in j:n
    ij += 1
    L[ij] = i+j
end
```
which loops over all lower triangle entries of `L::LowerTriangularArray` and the single
index `ij` is simply counted up. However, one could also use `[i, j]` as indices in the
loop body or to perform any calculation (`i+j` here).

!!! warn "`end` doesn't work for matrix indexing"
    Indexing `LowerTriangularMatrix` and `LowerTriangularArray` in matrix style (`[i, j]`) with 
    `end` doesn't work. It either returns an error or wrong results as the `end` is lowered by 
    Julia to the size of the underlying flat array dimension.

The `setindex!` functionality of matrixes will throw a `BoundsError` when trying to write
into the upper triangle of a `LowerTriangularArray`, for example
```@repl LowerTriangularArrays
L[2, 1] = 0    # valid index

L[1, 2] = 0    # invalid index in the upper triangle
```

But reading from it will just return a zero

```@repl LowerTriangularArrays
L[2, 3]     # in the upper triangle
```

Higher dimensional `LowerTriangularArray` can be indexed with multidimensional array indices 
like most other arrays types. Both the vector index and the matrix index for the lower 
triangle work as shown here
```@repl LowerTriangularArrays
L = rand(LowerTriangularArray{Float32}, 3, 3, 5)

L[2, 1] # second lower triangle element of the first lower triangle matrix 

L[2, 1, 1] # (2,1) element of the first lower triangle matrix 
```
The `setindex!` functionality follows accordingly. 

### Iterators

An iterator over all entries in the array can be created with `eachindex`
```@repl LowerTriangularArrays
L = rand(LowerTriangularArray, 5, 5, 5)
for ij in eachindex(L)
    # do something
end

eachindex(L)
```

In order to only loop over the harmonics (essentially the horizontal, ignoring other dimensions)
use `eachharmonic`
```@repl LowerTriangularArrays
eachharmonic(L)
```

If you only want to loop over the other dimensions use `eachmatrix`

```@repl LowerTriangularArrays
eachmatrix(L)
```

together they can be used as

```@repl LowerTriangularArrays
for k in eachmatrix(L)
    for lm in eachharmonic(L)
        L[lm, k]
    end
end
```

Note that `k` is a `CartesianIndex` that will loop over *all* other dimensions, whether there's only 1
(representing a 3D variable) or 5 (representing a 6D variable with the first two dimensions being a
lower triangular matrix).

## Linear algebra with `LowerTriangularArray`

The [LowerTriangularArrays](@ref lowertriangularmatrices) module's main purpose is not linear algebra,
and typical matrix operations will not work with `LowerTriangularMatrix` because it's treated as a vector
not as a matrix, meaning that the following will not work as expected

```@repl LowerTriangularArrays
L = rand(LowerTriangularMatrix{Float32}, 3, 3)
L * L
inv(L)
```

And many other operations that require `L` to be a `AbstractMatrix` which it isn't. In contrast, typical
vector operations like a scalar product between two "LowerTriangularMatrix" vectors does work


```@repl LowerTriangularArrays
L' * L
```

Summation with `sum` follows the flat, single index logic
```@repl 
L = rand(LowerTriangularMatrix{Float32}, 3, 3, 5)
sum(L, dims=2) 
```
sums along the second dimension of the underlying vector, not of the full matrix representation. 

## Rotation of `LowerTriangularArray`

`LowerTriangularArray`s are used to describe spherical harmonics. In that case each element
represents the coefficient in fron of the respective harmonic describing a field on the sphere
when transformed to grid space. We implement `rotate!` (and `rotate` for an allocating version)
for `LowerTriangularArray` to rotate these coefficients in complex number space to represent
a longitude rotation of the represented grid space field. For example start with

```@repl LowerTriangularArrays
M = rand(LowerTriangularMatrix{ComplexF32}, 3, 3)
```

Now `rotate!(::LowerTriangularArray, degree)`

```@repl LowerTriangularArrays
rotate!(M, 45)
```

represents the same (up to rounding errors from the rotation when not rotating by ``\pm 90, \pm 180, ...``)
field in grid space but rotated by 45˚ eastward. Note how the zonal modes (the first column) aren't
rotated because they are zonally constant anyway (in fact their imaginary part can be dropped,
but isn't here as created by the `rand`) and for the other modes this amounts to a multiplication with

```math
\exp(-i\frac{2π}{360}\phi m)
```
With ``\phi`` the rotation angle in degrees (positive eastward) and ``m`` the zonal wavenumber (order of the
spherical harmonic, the zero-based column index). Rotating again by 315˚ yields the original array

```@repl LowerTriangularArrays
rotate!(M, 315)
```

## Reverse of `LowerTriangularArray`

A `LowerTriangularArray` is an `AbstractArray`, as such `reverse` and `reverse!` (in-place) are defined,
reversing all elements of these arrays in the way how they are indexed with a single index.
For `LowerTriangularArra` representing the coefficients of the spherical harmonics this does not make
much sense, however, we describe here the functionality to `reverse` these arrays as they represent
spherical harmonics, adding methods for `dims=:lat` for reversal in latitude direction and `dims=:lon` in
longitude direction. Spherical harmonics are reversed in latitude by flipping the sign of the
odd harmonics, which are the ones that are not symmetric around the equator. Spherical harmonics are
reversed in longitude by taking the complex conjugate of every element as this flips the sign of the
imaginary parts, which effectively mirrors the rotation of that harmonic around 0˚E.

```@repl LowerTriangularArrays
reverse(M, dims=:lat)
```

and in longitude

```@repl LowerTriangularArrays
reverse(M, dims=:lon)
```

## Broadcasting with `LowerTriangularArray`

In contrast to linear algebra, many element-wise operations work as expected thanks to broadcasting,
so operations that can be written in `.` notation whether implicit `+`, `2*`, ... or explicitly
written `.+`, `.^`, ... or via the `@.` macro

```@repl LowerTriangularArrays
L + L
2L

@. L + 2L - 1.1*L / L^2
```

## GPU 

`LowerTriangularArray{T, N, S, ArrayType}` wraps around an array of type `ArrayType`.
If this array is a GPU array (e.g. `CuArray`), all operations are performed on GPU as well (work in progress).
The implementation was written so that scalar indexing is avoided in almost all cases,
so that GPU operation should be performant.
To use `LowerTriangularArray` on GPU you can e.g. just `adapt` an existing `LowerTriangularArray`.

```julia 
using Adapt
L = rand(LowerTriangularArray{Float32}, 5, 5, 5)
L_gpu = adapt(CuArray, L)
```

## Function and type index

```@autodocs
Modules = [SpeedyWeather.LowerTriangularArrays]
```
