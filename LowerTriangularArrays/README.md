# LowerTriangularArrays.jl

[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/stable/lowertriangularmatrices/)
[![docs](https://img.shields.io/badge/documentation-main-blue.svg)](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/lowertriangularmatrices/)

LowerTriangularArrays is a package that has been developed for SpeedyWeather.jl but can also be used standalone. It's seperately registered as well, so you can install it via easily `using Pkg; Pkg.add("LowerTriangularArrays")`.

This module defines `LowerTriangularArray`, a lower triangular matrix format, which in contrast to
[`LinearAlgebra.LowerTriangular`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.LowerTriangular)
does not store the entries above the diagonal.
SpeedyWeather.jl uses `LowerTriangularArray` which is defined as a subtype of `AbstractArray` to store
the spherical harmonic coefficients. 
For 2D `LowerTriangularArray` the alias `LowerTriangularMatrix` exists.
Higher dimensional `LowerTriangularArray` are 'batches' of 2D `LowerTriangularMatrix`.
So, for example a ``(10\times 10\times 10)`` `LowerTriangularArray` holds 10 `LowerTriangularMatrix` of size ``(10\times 10)`` in one array.  

## Example Use 

We can create a random or zero `LowerTriangularArray` as follows
```julia
using LowerTriangularArrays

L = rand(LowerTriangularArray{Float32}, 5, 5)
L_zero = zeros(LowerTriangularArray{Float32}, 5, 5)
```
this allocates a `LowerTriangularArray` which corresponds to a matrix of size ``5\times 5`` but only saves the actual 15 entries of the lower triangle. The allocated arrays can be indexed both with a matrix index `[l, m]` denoting column and row or a flat index `[lm]` denoting the position in the flattened lower triangle. So for example, with Julia saving arrays in column-major order, it follows that for the arrays allocated above

```julia 
L[2, 1] == L[2]
L[2, 2] == L[6]
```

holds. Both styles of indexing can be used, but 
in performance-critical code a single index should be used, as this directly maps
to the index of the underlying data vector. The matrix index is somewhat slower
as it first has to be converted to the corresponding single index. 

For more details on the usage of `LowerTriangularArray` see the [documentation](https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/lowertriangularmatrices/).

## Related Modules

LowerTriangularArrays.jl is designed to be used in conjuction with the other modules of the SpeedyWeather.jl ecosystem: 
 
- **RingGrids** - Spherical grid definitions and operations
- **SpeedyTransforms** - Spherical harmonic transforms
- **SpeedyWeather** - The full atmospheric modeling library
