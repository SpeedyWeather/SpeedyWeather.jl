# Development notes

## Set up

Clone SpeedyWeather.jl from GitHub and change the current working directory:

```
git clone https://github.com/milankl/SpeedyWeather.jl.git
cd SpeedyWeather.jl
```

To install Julia, download the latest version from the [Julia release page](https://julialang.org/downloads/) or, alternately, for easier cross-platform support and environment management, install through [Miniconda](https://docs.conda.io/en/latest/miniconda.html) with:

```
conda env create -f environment.yml
```

Then, activate the environment with `conda activate SpeedyWeather`.


## Install

To install SpeedyWeather.jl as a development package and [Revise.jl](https://github.com/timholy/Revise.jl) that automatically updates function definitions in a running Julia session, run:

```
julia -e 'import Pkg; Pkg.develop(path="."); Pkg.add("Revise")'
```


## Develop

Make sure to load Revise before importing anything else. For example:

```julia
using Revise, SpeedyWeather
run_speedy(n_days=30, trunc=63, Grid=OctahedralGaussianGrid, model=:shallowwater, output=true)
```  


## Versioning

This project uses [semantic versioning](https://semver.org/).
