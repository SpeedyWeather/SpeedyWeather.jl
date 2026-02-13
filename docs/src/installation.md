# Installation

SpeedyWeather.jl is registered in the Julia Registry. In most cases just open the
[Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and type
```julia
julia> using Pkg
julia> Pkg.add("SpeedyWeather")
```
or, equivalently, (`]` opens the package manager)
```julia
julia> ] add SpeedyWeather
```
which will automatically install the [latest release](https://github.com/SpeedyWeather/SpeedyWeather.jl/releases)
and all necessary dependencies. If you run into any troubles please raise an
[issue](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/new).

## Installing main

On the branch `main` we develop SpeedyWeather. This branch contains the latest
bug fixes and features that have not been released yet but there might also
be issues with `main` we have not fixed yet. If you need to install that latest
version from main do

```julia
julia> Pkg.add(url="https://github.com/SpeedyWeather/SpeedyWeather.jl", rev="main", subdir="SpeedyWeather")
```
`rev` (revision) refers to the `main` branch, and `subdir` is needed as we have structured
SpeedyWeather as a monorepo with other packages living in the same repository
(do `subdir="RingGrids"` for example if you want to install the latest unreleased version of
the RingGrids package ...). In a similar manner, you can also install other branches than `main`,
e.g. from a specific pull request.

In brief you can do the same as

```julia
(@v1.12) pkg> add SpeedyWeather:SpeedyWeather#main
```

following a `Repository:Subdirectory#branch` logic. Note that installing `main` gives you the
latest unreleased version of SpeedyWeather but the latest released versions of its subpackages.
When we move things between these packages then you may run into issues here. In that
case you need to `develop` the package as outlined in the next section.

## Developing SpeedyWeather

Julia's package manager Pkg.jl allows another way to develop/use the latest version
on main

```julia
julia> Pkg.develop(url="https://github.com/SpeedyWeather/SpeedyWeather.jl", subdir="SpeedyWeather")
```

in brief you can do this in the pkg shell with

```julia
(@v1.12) pkg> dev SpeedyWeather:SpeedyWeather
```

The first "SpeedyWeather" refers to the repository the second to the subdirectory
(which could be RingGrids, etc. as well) as also declared above.
See the [Managing packages](https://pkgdocs.julialang.org/dev/managing-packages/#Managing-Packages)
section in the Pkg documentation to learn more about the differences between `add` and `dev`.

## Compatibility with Julia versions

SpeedyWeather.jl requires Julia v1.10 or later. The package is tested on Julia 1.10, and 1.11.

## Extensions

SpeedyWeather.jl has a weak dependency on

- [Makie.jl](https://github.com/MakieOrg/Makie.jl) to extend `Makie.heatmap`

This is an
[extension](https://pkgdocs.julialang.org/v1.10/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)),
meaning that this functionality is only loaded from SpeedyWeather when `using Makie`
(or its backends CairoMakie.jl, GLMakie.jl, ...). Hence, if you want to make use of this
extension (some [Examples](@ref Examples) show this) you need to install Makie.jl manually.

