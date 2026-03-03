# Input data

Several components will allow you to specify your own input data. You
can either do this with the [`set!`](@ref) method, passing on a
`Field` to a variable, or many components allow for a `path` and `filename`.
See for example [Load orography from file](@ref).

However, we also use SpeedyWeatherAssets.jl to handle different versions
of input data with Julia's Artifacts system (`Artifacts.jl` and `Pkg.Artifacts`)
which will automatically download such data on-demand. This helps us to
reduce SpeedyWeather's repository size to a minimum while only downloading
data that is actually needed, and only downloading it once. For example,
this means that if you only use SpeedyWeather's `BarotropicModel` you
never have to download orography or a land-sea mask.

The following explains how we set up
[SpeedyWeatherAssets](https://github.com/SpeedyWeather/SpeedyWeatherAssets)
to work with `SpeedyWeather.get_asset` using different versions of
input data or from experimental branches. The primary intended use case
is to support different data-driven parameterizations with specific
neural network weights that all come in different files that need to be
loaded in.

## SpeedyWeatherAssets

The interface to load data from
[SpeedyWeatherAssets](https://github.com/SpeedyWeather/SpeedyWeatherAssets)
is via the `get_asset` function

```@docs; canonical=false
get_asset
```

with `from_assets = true` will look into that repository to find the
corresponding file under `path`, note that the base path starts with
`"data"` and then can have any subfolder(s). 

```julia
joinpath("data", "subfolder", file)
```

with `from_assets = false` will look for this `file` locally instead.

The function `get_asset` has a keyword argument `version`, use

- `v"1"` (a `VersionNumber`) for released versions of SpeedyWeatherAssets, we follow semantic versioning.
- `"branch"` (a `String`) to use data in a given `branch` of the SpeedyWeatherAssets repository.

In the latter case the version number (here `1`) is actually ignored and the

Then there are options how to read the data from file and interpret it.
`FileFormat` determines the file format, e.g. `NCDataset` will which attempt to
read in a netCDF file. `name` will read a variable of that name inside the file
and `ArrayType` is the type the data array is converted to or wrapped in, in most
cases this would be a `FullGaussianField` or `FullClenshawField`.
Note that this excludes the `on_architecture` call that will automatically convert
to a `GPUArray`. To define a new way to read in and interpret an asset extend
the following function

```julia
_get_asset(path::String, name::String, ArrayType::Type{<:YourArraType}, FileFormat::Type{<:YourFileFormat})
```

for `YourArrayType` and `YourFileFormat`.

When calling `get_asset` a new asset will be downloaded and stored where
other Julia Artifacts are stored, by default this would be `.julia/artifacts`
and then folders of their respective hash. Using assets in SpeedyWeather one
does not handle the Artifacts.toml file directly but this will be placed in the
root folder (or inside `SpeedyWeather/src/input` when loading locally with `.SpeedyWeather`).
Deleting this file will redownload all demanded assets but one can also
delete individiual ones. 

## Managing new SpeedyWeatherAssets versions

The current default version of SpeedyWeatherAssets used is

```@example assets
using SpeedyWeather
SpeedyWeather.DEFAULT_ASSETS_VERSION
```

For new released versions one should (manually) add that version to

```@example assets
SpeedyWeather.AVAILABLE_ASSETS_VERSIONS
```

only those versions can be used in `version = v"1.1"` for example.
Branches of SpeedyWeatherAssets, e.g. `version = "main"` do not have
to be added to the available versions.



