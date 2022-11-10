# How to run SpeedyWeather.jl

The simplest way to run SpeedyWeather.jl with default parameters is

```julia
using SpeedyWeather
run_speedy()
```

Hooray, you have just simulated the Earth's atmosphere. Parameters, their meanings and
defaults can be found in `src/default_parameters.jl`, an incomplete list is provided below.
For example, if you want the simulation to run in double precision (`Float64`), at higher
resolution (`trunc`, the triangular spectral truncation), slow down the rotation of the Earth
(`rotation_earth` in ``s^{-1}``) and create some netCDF ouput, do

```julia
run_speedy(Float64,trunc=85,rotation_earth=1e-5,output=true)
```

If provided, the number format has to be the first argument, all other arguments are keyword arguments.

## Incomplete list of parameters

The following list only contains some frequently used parameters. As SpeedyWeather.jl is
developed this is by no means complete.


### Resolution

| parameter keyword | default value | meaning   |
| ----------------- | ------------- | --------- |
| `NF::DataType`    | `Float32`     | Number format (not keyword, first positional, but optional argument) |
| `trunc::Int`      | `31`          | Spectral triangular truncation , i.e. the maximum order and degree of the spherical harmonics | 
| `nlev::Int`       | `8`           | Number of vertical levels |


### Physical parameters

| parameter keyword | default value | meaning   |
| ----------------- | ------------- | --------- |
| `model::Type{<:ModelSetup}`   | `BarotropicModel` | Equations to be solved: `BarotropicModel`, `ShallowWaterModel`, or `PrimitiveEquationModel` |
| `radius_earth::Real`  | `6.371e6` | Radius of Earth [``m``] |
| `rotation_earth::Real`| `7.292e-5`| Angular frequency of Earth's rotation [``s^{-1}``] |
| `gravity::Real`   | `9.81`        | Gravitational acceleration [``m/s^2``] |


### Diffusion

| parameter keyword | default value | meaning   |
| ----------------- | ------------- | --------- |
| `diffusion_power::Real` | `4`     | Power `n` of Laplacian in horizontal diffusion ``\nabla^{2n}`` |
| `diffusion_time::Real` | `2.4`    | Diffusion time scale [``hrs``] for temperature and vorticity |
| `diffusion_time_div::Real` | `diffusion_time` | Diffusion time scale [``hrs``] for divergence |      
| `diffusion_time_strat::Real` | `12`           | Diffusion time scale [``hrs``] for extra ``\nabla^2`` in the stratosphere |
| `damping_time_strat::Real` | `720`            | Damping time [``hrs``] for drag on zonal-mean wind in the stratosphere |


### Time stepping

| parameter keyword | default value | meaning   |
| ----------------- | ------------- | --------- |
| `Δt_at_T31::Real` | `60`          | Time step in minutes for T31, scale linearly for specified `trunc` |
| `n_days::Real`    | `10`          | Number of days to integrate for |
| `robert_filter::Real` | `0.05`    | Robert (1966) time filter coefficeint to suppress computational mode |
| `williams_filter::Real` | `0.53`  | William's time filter (Amezcua 2011) coefficient for 3rd order acc |


### Spectral transforms

| parameter keyword | default value | meaning   |
| ----------------- | ------------- | --------- |
| `recompute_legendre::Bool` | `false` | Recompute or precompute the Legendre polynomials in the transforms |


### Output

| parameter keyword | default value | meaning   |
| ----------------- | ------------- | --------- |
| `verbose::Bool`   | `true`        | Print dialog for feedback |
| `output::Bool`    | `false`       | Store data in netCDF? |
| `output_dt::Real` | `6`           | Output time step [``hrs``] |
| `startdate::DateTime` |  `DateTime(2000,1,1)` | Start date of the time axis (default: Jan 1, 2000) |
| `output_path::String`| `pwd()`       | Path to output folder (default is current directory) |


## The `run_speedy` interface

```@docs
run_speedy
```

## The `initialize_speedy` interface

```@docs
initialize_speedy
```