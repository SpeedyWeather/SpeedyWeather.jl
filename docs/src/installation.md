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

However, you may want to make use of the latest features, then install directly from the main branch with
```julia
julia> Pkg.add(url="https://github.com/SpeedyWeather/SpeedyWeather.jl", rev="main")
```
In a similar manner, other branches can be also installed. You can also type, equivalently,
```julia
julia> ] add https://github.com/SpeedyWeather/SpeedyWeather.jl#main
```

## Compatibility with Julia versions

SpeedyWeather.jl requires Julia v1.9 or later. The package is tested on Julia 1.9, and 1.10.

## Extensions

SpeedyWeather.jl has a weak dependency on

- [Makie.jl](https://github.com/MakieOrg/Makie.jl) to extend `Makie.heatmap`

This is an
[extension](https://pkgdocs.julialang.org/v1.10/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)),
meaning that this functionality is only loaded from SpeedyWeather when `using Makie`
(or its backends CairoMakie.jl, GLMakie.jl, ...). Hence, if you want to make use of this
extension (some [Examples](@ref Examples) show this) you need to install Makie.jl manually.

