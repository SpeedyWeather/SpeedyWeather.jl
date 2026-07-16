"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ZonalVelocity10mOutput <: AbstractOutputVariable
    name::String = "u10"
    unit::String = "m/s"
    long_name::String = "10m zonal wind speed"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

path(::ZonalVelocity10mOutput, simulation) = simulation.variables.grid.u

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MeridionalVelocity10mOutput <: AbstractOutputVariable
    name::String = "v10"
    unit::String = "m/s"
    long_name::String = "10m meridional wind speed"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

path(::MeridionalVelocity10mOutput, simulation) = simulation.variables.grid.v

"""$(TYPEDSIGNATURES)
10m wind is defined using a logarithmic profile from the lowermost model layer.

    u_10 = u_bottom * ln(10 / z₀) / ln(z_bottom / z₀)

(with `z` in m), with

    z_bottom = z_surf + T_bottom * Δp_geopot_bottom / g
"""
function output!(
        output::AbstractOutput,
        variable::Union{ZonalVelocity10mOutput, MeridionalVelocity10mOutput},
        simulation::AbstractSimulation,
    )
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    # Retrieve u_bottom
    var = path_or_nothing(variable, simulation)
    isnothing(var) && return nothing                # silently escape early if variable is not defined in the simulation
    u_or_v_grid = get_prognostic_step(var, simulation.model.time_stepping, output)
    nlayers = size(u_or_v_grid, 2)
    u_or_v_bottom = field_view(u_or_v_grid, :, nlayers)

    z_bottom = simulation.variables.scratch.grid.a_2D
    u_or_v10 = simulation.variables.scratch.grid.b_2D

    # Compute z_bottom as z_surf + T_bottom * R/g * log(pₛ/p_full_bottom), start with z_surf
    z_bottom .= max.(simulation.model.orography.orography, 0)   # [m] set negative values to zero
    T = get_prognostic_step(simulation.variables.grid.temperature, simulation.model.time_stepping, output)
    T_bottom = field_view(T, :, nlayers)
    R_dry = simulation.model.atmosphere.R_dry
    g = simulation.model.planet.gravity
    coord = simulation.model.geometry.vertical_coordinates
    pₛ = simulation.variables.parameterizations.surface_pressure

    # R/g * log(pₛ/p_full_bottom): for sigma coords this is constant (= Δp_geopot_full[end]/g);
    # for hybrid it varies per grid point since p_full_bottom = A*p_ref + B*pₛ.
    @. z_bottom += T_bottom * (R_dry / g) * log(pₛ / pressure(nlayers, pₛ, coord))

    # Compute u10, TODO should this be the same z₀ as in vertical diffusion or surface fluxes?
    z₀ = simulation.variables.parameterizations.surface_roughness

    # include a max here as this will throw an error if the model blows up and have negative temperatures
    # if the max case is hit we divide by zero which yields inf so clearly flagged in any case
    @. u_or_v10 = u_or_v_bottom .* log(10 / z₀) ./ log.(max.(z_bottom, z₀) / z₀)

    # interpolate 2D/3D variables
    u_or_v10_output = output.field2D
    u_or_v10_grid = on_architecture(CPU(), u_or_v10)
    RingGrids.interpolate!(u_or_v10_output, u_or_v10_grid, output.interpolator)

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(u_or_v10_output, variable.keepbits)
    end

    write_array!(output, variable, u_or_v10_output)
    return nothing
end

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceTemperatureOutput{F} <: AbstractOutputVariable
    name::String = "tsurf"
    unit::String = "degC"
    long_name::String = "Surface air temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
    transform::F = (x) -> x - 273.15     # [K] to [˚C]
end

# not the actual surface temperature but the core variable to read in
path(::SurfaceTemperatureOutput, simulation) = simulation.variables.grid.temperature

function output!(
        output::AbstractOutput,
        variable::SurfaceTemperatureOutput,
        simulation::AbstractSimulation,
    )
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    # reuse scratch array to avoid allocations
    Ts = simulation.variables.scratch.grid.a_2D
    pN = simulation.variables.scratch.grid.b_2D

    # Retrieve T_bottom
    var = path_or_nothing(variable, simulation)
    isnothing(var) && return nothing                # silently escape early if variable is not defined in
    T_grid = get_prognostic_step(var, simulation.model.time_stepping, output)
    nlayers = size(T_grid, 2)
    T_bottom = field_view(T_grid, :, nlayers)

    # Other parameters
    κ = simulation.model.atmosphere.κ
    coord = simulation.model.geometry.vertical_coordinates
    pₛ = simulation.variables.parameterizations.surface_pressure
    pN .= pressure.(nlayers, pₛ, coord)

    # Compute Ts assuming dry adiabatic profile: T_surf = T_bottom * (pₛ / p_bottom)^κ
    (; transform) = variable
    @. Ts = transform(T_bottom * (pₛ / pN)^κ)   # Convert to °C

    # interpolate 2D/3D variables
    Ts_output = output.field2D
    Ts_grid = on_architecture(CPU(), Ts)
    RingGrids.interpolate!(Ts_output, Ts_grid, output.interpolator)

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(Ts_output, variable.keepbits)
    end

    write_array!(output, variable, Ts_output)
    return nothing
end


"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MomentumBoundaryLayerDragOutput <: AbstractOutputVariable
    name::String = "bld_m"
    unit::String = "1"
    long_name::String = "Momentum flux related boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::MomentumBoundaryLayerDragOutput, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag_momentum

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct HeatBoundaryLayerDragOutput <: AbstractOutputVariable
    name::String = "bld_h"
    unit::String = "1"
    long_name::String = "Heat flux related boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::HeatBoundaryLayerDragOutput, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag_heat

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MoistureBoundaryLayerDragOutput <: AbstractOutputVariable
    name::String = "bld_q"
    unit::String = "1"
    long_name::String = "Moisture flux related boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::MoistureBoundaryLayerDragOutput, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag_moisture

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MomentumSurfaceRoughnessOutput <: AbstractOutputVariable
    name::String = "sr"
    unit::String = "m"
    long_name::String = "Surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::MomentumSurfaceRoughnessOutput, simulation) =
    simulation.variables.parameterizations.momentum_roughness

@kwdef mutable struct HeatSurfaceRoughnessOutput <: AbstractOutputVariable
    name::String = "sr_h"
    unit::String = "m"
    long_name::String = "Surface heat roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::HeatSurfaceRoughnessOutput, simulation) =
    simulation.variables.parameterizations.heat_roughness

# collect all in one for convenience
BoundaryLayerOutput() = (
    ZonalVelocity10mOutput(),
    MeridionalVelocity10mOutput(),
    SurfaceTemperatureOutput(),
    MomentumBoundaryLayerDragOutput(),
    HeatBoundaryLayerDragOutput(),
    MoistureBoundaryLayerDragOutput(),
    MomentumSurfaceRoughnessOutput(),
    HeatSurfaceRoughnessOutput(),
)
