# [Array dimensions](@id array_dimensions)

SpeedyWeather's two main array types, `Field` (see [RingGrids](@ref)) for gridded data
and `LowerTriangularArray` (see [LowerTriangularArrays](@ref lowertriangularmatrices))
for spectral coefficients, unravel the horizontal into the first array dimension.
Any additional array dimension can then represent whatever you like: the vertical,
time (or time steps of the time integration), etc. To record what these additional
dimensions actually mean, both types (optionally) carry a dimension tag from the
`ArrayDimensions` module alongside their data. The following tags are defined

|                                   | 2D (horizontal only) | 3D + vertical | 3D + time | 4D + vertical + time |
|-----------------------------------|----------------------|---------------|-----------|----------------------|
| grid (`Field`)                    | `XY` (default)       | `XYZ`         | `XYT`     | `XYZT`               |
| spectral (`LowerTriangularArray`) | `LM` (default)       | `LMZ`         | `LMT`     | `LMZT`               |

with `X`, `Y` for longitude, latitude and `L`, `M` for degree, order of the spherical
harmonics (both pairs unravelled into a single array dimension), `Z` for the vertical and
`T` for time. For a `ColumnField` (vertical first) there are also `ZXY` and `ZXYT`.
The defaults `XY` and `LM` make no assumption about the meaning of any additional dimension.

All of these are singleton types: they hold no data, do not change the memory layout of
the arrays, and, being known at compile time, come at no runtime cost.

## Where to see the dimensions

The array summary includes the dimension tag in parentheses

```@example dimensions
using SpeedyWeather
spectrum = Spectrum(trunc=5)
L = rand(ComplexF32, spectrum, ArrayDimensions.LMZ(), 3)
```

here `(LMZ)` says that the 3 elements in the second array dimension are vertical layers.
Without an explicit tag the default applies

```@example dimensions
L2 = rand(ComplexF32, spectrum, 3)
```

and the summary reads `(LM+)`: the `+` denotes that the array has more dimensions than
the tag describes, i.e. no assumption is made about what the second dimension means.
The tag is also directly accessible as a field

```@example dimensions
L.dims
```

All of the above works the same for `Field`, e.g.

```@example dimensions
field = zeros(FullGaussianGrid(4), ArrayDimensions.XYZ(), 3)
```

## What the dimensions are (not) used for

The dimension tags are bookkeeping only. They do not change how an array is indexed,
computed on, or broadcast: multiplying an `XYZ` field with an `XYT` field will *not*
return an `XYZT` field — broadcasting ignores the tags entirely (the result simply
carries over the tag of one of its inputs). Instead, the tags allow manual decisions
whenever the meaning of a dimension matters, through

```@example dimensions
ArrayDimensions.hasvertical(L), ArrayDimensions.hastime(L)
```

or through dispatch on the tag types and their unions like
`ArrayDimensions.DimensionsWithTime` and `ArrayDimensions.DimensionsWithVertical`.

The tags are preserved through `similar`, `zero`, views and broadcasting; indexing
into a tagged dimension with an integer drops it accordingly, e.g. `L[:, 1]` of an
`LMZ`-tagged array returns an `LM`-tagged one.

Within SpeedyWeather, the variables in `simulation.variables` carry these tags; in
particular the step dimension that prognostic and tendency variables have for the time
integration is tagged as time `T`, so their summaries show, e.g., `(LMZT)` or `(XYZT)`,
see [Step dimension](@ref). For more details on the two array types see
[Array dimensions of a Field](@ref) and
[Array dimensions for `LowerTriangularArray`](@ref).
