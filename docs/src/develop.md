# Development notes

## Set up

Clone SpeedyWeather.jl from GitHub and change the current working directory:

```
git clone https://github.com/milankl/SpeedyWeather.jl.git
cd SpeedyWeather.jl
```

To install Julia, manually download the latest version from the [Julia release page](https://julialang.org/downloads/) or, alternately, install it using [juliaup](https://github.com/JuliaLang/juliaup), the official Julia installer and version multiplexer. For Linux, macOS, and WSL run:

```
curl -fsSL https://install.julialang.org | sh
```


## Install

To install SpeedyWeather.jl as a development package and [Revise.jl](https://github.com/timholy/Revise.jl) to automatically update function definitions in a running session, run:

```
julia -e 'import Pkg; Pkg.develop(path="."); Pkg.add("Revise")'
```


## Develop

Make sure to load Revise before importing anything else. For example:

```julia
using Revise, SpeedyWeather
run_speedy(n_days=30, trunc=63, Grid=OctahedralGaussianGrid, model=:shallowwater, output=true)
```  

## Documentation

To generate the documentation, run:

```
julia --project=docs -e 'import Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```


## Testing

To run tests, run: 

```
julia -e 'import Pkg; Pkg.test("SpeedyWeather")'
```


## Versioning

This project uses [semantic versioning](https://semver.org/).
