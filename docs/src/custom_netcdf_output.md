# Customizing netCDF output

SpeedyWeather's [NetCDF output](@ref) is modularised for the output variables,
meaning you can add relatively easy new variables to be outputted
alongside the default variables in the netCDF file. We explain here
how to define a new output variable largely following the logic
of [Extending SpeedyWeather](@ref).

## Define a new output variable

Say we want to output the [Vertical velocity](@ref). In [Sigma coordinates](@ref)
on every time step, one has to integrate the divergence vertically to
know where the flow is not divergence free, meaning that the horizontally
converging or diverging motion is balanced by a vertical velocity.
This leads to the variable ``\partial \sigma / \partial t``, which
is the equivalent of [Vertical velocity](@ref) in the [Sigma coordinates](@ref).
This variable is calculated and stored at every time step in 

```julia
simulation.diagnostic_variables.dynamics.σ_tend
```

So how do we access it and add it the netCDF output?

First we define `VerticalVelocityOutput` as a new `struct` subtype of
`SpeedyWeather.AbstractOutputVariable` we add the required fields
`name::String`, `unit::String`, `long_name::String` and
`dims_xyzt::NTuple{4, Bool}` (we skip the optional fields
for `missing_value` or compression).

```@example netcdf_custom
using SpeedyWeather

@kwdef struct VerticalVelocityOutput <: SpeedyWeather.AbstractOutputVariable
    name::String = "w"
    unit::String = "s^-1"
    long_name::String = "vertical velocity"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
end
```

By default (using the `@kwdef` macro) we set the dimensions in `dims_xyzt`
to 4D because the vertical velocity is a 3D variable that we want to output
on every time step. So while it is a required field for every output variable
you should not actually change it as it is an inherent property of the output
variable.

You can now add this variable to the `NetCDFOutput` as already described in
[Output variables](@ref)

```@example netcdf_custom
spectral_grid = SpectralGrid()
output = NetCDFOutput(spectral_grid)
add!(output, VerticalVelocityOutput())
```

Note that here we skip the `SpeedyWeather.` which would point to the
SpeedyWeather scope as we have defined `VerticalVelocityOutput` in
the global scope.

## Define how to output a new variable

While we have defined a new output variable we have not actually
defined how to output it. Because in the end we will need to
write that variable into the netcdf file in `NetCDFOutput`,
which we describe now.

```@example netcdf_custom
function output!(
    output::NetCDFOutput,
    variable::VerticalVelocityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    w = output.grid3D
    (; σ_tend) = diagn.dynamics
    RingGrids.interpolate!(w, σ_tend , output.interpolator)

    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = w
    return nothing
end
```