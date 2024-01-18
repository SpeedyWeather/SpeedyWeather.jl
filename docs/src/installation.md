# Installation

SpeedyWeather.jl is registered in the Julia Registry. In most cases just open the
[Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and type
```julia
julia> using Pkg
julia> Pkg.add("SpeedyWeather")
```
or, equivalently, (`]` opens the package manager)
```julia
julia>] add SpeedyWeather
```
which will automatically install the [latest release](https://github.com/SpeedyWeather/SpeedyWeather.jl/releases)
and all necessary dependencies. If you run into any troubles please raise an
[issue](https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/new).

However, you may want to make use of the latest features, then install directly from the main branch with
```julia
julia> Pkg.add(url="https://github.com/SpeedyWeather/SpeedyWeather.jl",rev="main")
```
other branches than `main` can be similarly installed. You can also type, equivalently,
```julia
julia>] add https://github.com/SpeedyWeather/SpeedyWeather.jl#main
```

## Compatibility with Julia versions

SpeedyWeather.jl requires Julia v1.8 or later. The package is tested on Julia v1.8, 1.9, and 1.10.
