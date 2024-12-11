## VORTICITY -------------

"""Defines netCDF output of vorticity. Fields are
$(TYPEDFIELDS)

Custom variable output defined similarly with required fields marked,
optional fields otherwise use variable-independent defaults. Initialize with `VorticityOutput()`
and non-default fields can always be passed on as keyword arguments,
e.g. `VorticityOutput(long_name="relative vorticity", compression_level=0)`."""
@kwdef mutable struct VorticityOutput <: AbstractOutputVariable

    "[Required] short name of variable (unique) used in netCDF file and key for dictionary"
    name::String = "vor"

    "[Required] unit of variable"
    unit::String = "s^-1"

    "[Required] long name of variable used in netCDF file"
    long_name::String = "relative vorticity"

    "[Required] NetCDF dimensions the variable uses, lon, lat, layer, time"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)

    "[Optional] missing value for the variable, if not specified uses NaN"
    missing_value::Float64 = NaN

    "[Optional] compression level of the lossless compressor, 1=lowest/fastest, 9=highest/slowest, 3=default"
    compression_level::Int = 3

    "[Optional] bitshuffle the data for compression, false = default"
    shuffle::Bool = true

    "[Optional] number of mantissa bits to keep for compression (default: 15)"
    keepbits::Int = 5
end

"""$(TYPEDSIGNATURES)
Output the vorticity field `vor` from `diagn.grid` into the netCDF file `output.netcdf_file`.
Interpolates the vorticity field onto the output grid and resolution as specified in `output`.
Method required for all output variables `<: AbstractOutputVariable` with dispatch over the
second argument."""
function output!(
    output::NetCDFOutput,
    variable::VorticityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    vor = output.grid3D             # use output grid, vorticity is 3D
    (; vor_grid) = diagn.grid       # unpack vorticity from gridded diagnostic variables
    RingGrids.interpolate!(vor, vor_grid, output.interpolator)

    unscale!(vor, diagn.scale[])    # was vor*radius, back to vor
    round!(vor, variable.keepbits)  # bitrounding for compression
    i = output.output_counter       # write into timestep i
    output.netcdf_file[variable.name][:, :, :, i] = vor
    return nothing
end

## U velocity -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct ZonalVelocityOutput <: AbstractOutputVariable
    name::String = "u"
    unit::String = "m/s"
    long_name::String = "zonal wind"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

"""$(TYPEDSIGNATURES)
`output!` method for `ZonalVelocityOutput` to write the zonal velocity field `u` from `diagn.grid`,
see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::ZonalVelocityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    u = output.grid3D
    (; u_grid) = diagn.grid
    RingGrids.interpolate!(u, u_grid, output.interpolator)

    round!(u, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = u
    return nothing
end

## V velocity -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct MeridionalVelocityOutput <: AbstractOutputVariable
    name::String = "v"
    unit::String = "m/s"
    long_name::String = "meridional wind"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::MeridionalVelocityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    v = output.grid3D
    (; v_grid) = diagn.grid
    RingGrids.interpolate!(v, v_grid, output.interpolator)

    round!(v, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = v
    return nothing
end

## DIVERGENCE -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct DivergenceOutput <: AbstractOutputVariable
    name::String = "div"
    unit::String = "s^-1"
    long_name::String = "divergence"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 5
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::DivergenceOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    div = output.grid3D
    (; div_grid) = diagn.grid
    RingGrids.interpolate!(div, div_grid, output.interpolator)

    unscale!(div, diagn.scale[])    # was vor*radius, back to vor
    round!(div, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = div
    return nothing
end

## INTERFACE DISPLACEMENT -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct InterfaceDisplacementOutput <: AbstractOutputVariable
    name::String = "eta"
    unit::String = "m"
    long_name::String = "interface displacement"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::InterfaceDisplacementOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    eta = output.grid2D
    (; pres_grid) = diagn.grid
    RingGrids.interpolate!(eta, pres_grid, output.interpolator)

    round!(eta, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = eta
    return nothing
end

## SURFACE PRESSURE -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct SurfacePressureOutput <: AbstractOutputVariable
    name::String = "pres"
    unit::String = "hPa"
    long_name::String = "surface pressure"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::SurfacePressureOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    pres = output.grid2D
    (; pres_grid) = diagn.grid
    RingGrids.interpolate!(pres, pres_grid, output.interpolator)

    @inbounds for ij in eachindex(pres)
        pres[ij] = exp(pres[ij]) / 100    # from log(Pa) to hPa
    end

    round!(pres, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, i] = pres
    return nothing
end

## TEMPERATURE -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct TemperatureOutput <: AbstractOutputVariable
    name::String = "temp"
    unit::String = "degC"
    long_name::String = "temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::TemperatureOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    temp = output.grid3D
    (; temp_grid) = diagn.grid
    RingGrids.interpolate!(temp, temp_grid, output.interpolator)
    temp .-= 273.15             # convert from K to ˚C

    round!(temp, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = temp
    return nothing
end

## HUMIDITY -------------

"""Defines netCDF output for a specific variables, see `VorticityOutput` for details.
Fields are $(TYPEDFIELDS)"""
@kwdef mutable struct HumidityOutput <: AbstractOutputVariable
    name::String = "humid"
    unit::String = "kg/kg"
    long_name::String = "specific humidity"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

"""$(TYPEDSIGNATURES)
`output!` method for `variable`, see `output!(::NetCDFOutput, ::VorticityOutput, ...)` for details."""
function output!(
    output::NetCDFOutput,
    variable::HumidityOutput,
    progn::PrognosticVariables,
    diagn::DiagnosticVariables,
    model::AbstractModel,
)
    humid = output.grid3D
    (; humid_grid) = diagn.grid
    RingGrids.interpolate!(humid, humid_grid, output.interpolator)

    round!(humid, variable.keepbits)
    i = output.output_counter   # output time step to write
    output.netcdf_file[variable.name][:, :, :, i] = humid
    return nothing
end