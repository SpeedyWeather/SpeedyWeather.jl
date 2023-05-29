# Installation

SpeedyWeather.jl is registered in the Julia Registry. In most cases just open the
[Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and type
```julia
julia> using Pkg
julia> Pkg.add("SpeedyWeather")
```
which will automatically install the [latest release](https://github.com/SpeedyWeather/SpeedyWeather.jl/releases)
and all necessary dependencies. If you run into any troubles please raise an
[issue](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/new)

However, you may want to make use of the latest features, then install directly from the main branch with
```julia
julia> Pkg.add(url="https://github.com/SpeedyWeather/SpeedyWeather.jl",rev="main")
```
other branches than `main` can be similarly installed.

## Compatibility with Julia versions

SpeedyWeather.jl usually lives on the latest minor release and/or its predecessor.
At the moment (May 2023) this means 
- Julia v1.8
- Julia v1.9

are supported, but we dropped the support of earlier versions.
