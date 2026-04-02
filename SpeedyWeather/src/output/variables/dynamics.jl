"""Defines netCDF output of vorticity. Fields are
$(TYPEDFIELDS)

Custom variable output defined similarly with required fields marked,
optional fields otherwise use variable-independent defaults. Initialize with `VorticityOutput()`
and non-default fields can always be passed on as keyword arguments,
e.g. `VorticityOutput(long_name="relative vorticity", compression_level=0)`.
Custom variable output also requires the `path(::MyOutputVariable, simulation)`
to be extended to return the AbstractField subject to output.
Custom element-wise variable transforms, e.g. scale and/or offset to change
units, or even exp(x)/100 to change from log surface pressure to hPa
are passed on as `transform::Function = x -> exp(x)/100`."""
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

    "[Optional] Unscale the variable for output? (default: true)"
    unscale::Bool = true
end

"""$TYPEDSIGNATURES To be extended for every output variable to define
the path where in `simulation` to find that output variable `::AbstractField`."""
path(::VorticityOutput, simulation) = simulation.variables.grid.vorticity

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
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

path(::ZonalVelocityOutput, simulation) = simulation.variables.grid.u

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
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

path(::MeridionalVelocityOutput, simulation) = simulation.variables.grid.v

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct DivergenceOutput <: AbstractOutputVariable
    name::String = "div"
    unit::String = "s^-1"
    long_name::String = "divergence"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 5
    unscale::Bool = true
end

path(::DivergenceOutput, simulation) = simulation.variables.grid.divergence

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
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

path(::InterfaceDisplacementOutput, simulation) = simulation.variables.grid.η

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfacePressureOutput{F} <: AbstractOutputVariable
    name::String = "pres"
    unit::String = "hPa"
    long_name::String = "surface pressure"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
    transform::F = (x) -> exp(x) / 100     # log(Pa) to hPa
end

path(::SurfacePressureOutput, simulation) = simulation.variables.grid.pressure

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MeanSeaLevelPressureOutput{F} <: AbstractOutputVariable
    name::String = "mslp"
    unit::String = "hPa"
    long_name::String = "mean sea-level pressure"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
    transform::F = (x) -> exp(x) / 100     # log(Pa) to hPa
end

# points to surface not mean sea level pressure but core variable to read in
path(::MeanSeaLevelPressureOutput, simulation) = simulation.variables.grid.pressure

function output!(
        output::NetCDFOutput,
        variable::MeanSeaLevelPressureOutput,
        simulation::AbstractSimulation,
    )
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    # get log(surface pressure) field
    lnpₛ = path(variable, simulation)
    h = simulation.model.orography.orography
    (; R_dry, κ) = simulation.model.atmosphere
    g = simulation.model.planet.gravity

    # compute virtual temperature on the fly
    nlayers = size(simulation.variables.grid.temperature, 2)
    T = simulation.variables.parameterizations.surface_air_temperature

    if simulation.model.dynamics_only    # otherwise this has been computed already
        # calculate the surface air temperature from lowest model level temperature
        # via dry adiabatic lapse rate
        T .= field_view(simulation.variables.grid.temperature, :, nlayers)
        # σ vertical coordinate at lowest model level
        GPUArrays.@allowscalar σ = simulation.model.geometry.σ_levels_full[nlayers]
        σ⁻ᵏ = σ^(-κ)    # precalculate adiabatic descent factor
        T .*= σ⁻ᵏ       # lower to surface assuming dry adiabatic lapse rate
    end

    has_humid = haskey(simulation.variables.grid, :humidity)
    q = has_humid ? field_view(simulation.variables.grid.humidity, :, nlayers) : zero(T)
    Tᵥ = simulation.variables.scratch.grid.a_2D

    (; atmosphere) = simulation.model
    Tᵥ .= virtual_temperature.(T, q, atmosphere)

    # calculate mean sea-level pressure on model grid
    mslp = simulation.variables.scratch.grid.b_2D
    (; transform) = variable                    # to change units from log(Pa) to hPa
    @. mslp = transform(g * h / R_dry / Tᵥ + lnpₛ)    # Pa to hPa

    # interpolate 2D/3D variables
    mslp_output = output.field2D
    mslp_grid = on_architecture(CPU(), mslp)
    RingGrids.interpolate!(mslp_output, mslp_grid, output.interpolator)

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(mslp_output, variable.keepbits)
    end

    i = output.output_counter               # output time step i to write
    indices = get_indices(i, variable)      # returns (:, :, i) for example, depending on dims
    output.netcdf_file[variable.name][indices...] = mslp_output     # actually write to file
    return nothing
end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct TemperatureOutput{F} <: AbstractOutputVariable
    name::String = "temp"
    unit::String = "degC"
    long_name::String = "temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, true, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 10
    transform::F = (x) -> x - 273.15     # K to ˚C
end

path(::TemperatureOutput, simulation) = simulation.variables.grid.temperature

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
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

path(::HumidityOutput, simulation) = simulation.variables.grid.humidity

# collect all in one for convenience
DynamicsOutput() = (
    VorticityOutput(),
    ZonalVelocityOutput(),
    MeridionalVelocityOutput(),
    DivergenceOutput(),
    InterfaceDisplacementOutput(),
    SurfacePressureOutput(),
    MeanSeaLevelPressureOutput(),
    TemperatureOutput(),
    HumidityOutput(),
)
