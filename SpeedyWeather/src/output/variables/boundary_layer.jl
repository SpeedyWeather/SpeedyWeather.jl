"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct ZonalVelocity10m <: AbstractOutputVariable
    name::String = "u10"
    unit::String = "m/s"
    long_name::String = "10m zonal wind speed"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

path(::ZonalVelocity10m, simulation) = simulation.variables.grid.u

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MeridionalVelocity10m <: AbstractOutputVariable
    name::String = "v10"
    unit::String = "m/s"
    long_name::String = "10m meridional wind speed"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
end

path(::MeridionalVelocity10m, simulation) = simulation.variables.grid.v

"""$(TYPEDSIGNATURES)
10m wind is defined using a logarithmic profile from the lowermost model layer.

    u_10 = u_bottom * ln(10 / z₀) / ln(z_bottom / z₀)

(with `z` in m), with

    z_bottom = z_surf + T_bottom * Δp_geopot_bottom / g
"""
function interpolate_output!(
        field::AbstractField,
        output::AbstractOutput,
        variable::Union{ZonalVelocity10m, MeridionalVelocity10m},
        simulation::AbstractSimulation,
    )
    # Retrieve u_bottom
    u_or_v_grid = get_prognostic_step(path(variable, simulation), simulation.model.time_stepping, output)
    nlayers = size(u_or_v_grid, 2)
    u_or_v_bottom = field_view(u_or_v_grid, :, nlayers)

    z_bottom = simulation.variables.scratch.grid.a_2D
    u_or_v10 = simulation.variables.scratch.grid.b_2D

    # Compute z_bottom as z_surf + T_bottom * R/g * log(pₛ/p_full_bottom), start with z_surf
    z_bottom .= max.(simulation.model.orography.orography, zero(eltype(z_bottom)))   # [m] set negative values to zero
    T = get_prognostic_step(simulation.variables.grid.temperature, simulation.model.time_stepping, output)
    R_dry = simulation.model.atmosphere.R_dry
    g = simulation.model.planet.gravity
    coord = simulation.model.geometry.vertical_coordinates
    pₛ = simulation.variables.parameterizations.surface_pressure

    # R/g * log(pₛ/p_full_bottom): for sigma coords this is constant (= Δp_geopot_full[end]/g);
    # for hybrid it varies per grid point since p_full_bottom = A*p_ref + B*pₛ.
    # In a kernel so that the layer index stays an integer (Reactant traces numbers in broadcasts)
    launch!(
        architecture(z_bottom), LinearWorkOrder, size(z_bottom), bottom_layer_height_kernel!,
        z_bottom, T, pₛ, coord, nlayers, R_dry / g
    )

    # Compute u10, TODO should this be the same z₀ as in vertical diffusion or surface fluxes?
    z₀ = simulation.variables.parameterizations.momentum_roughness

    # include a max here as this will throw an error if the model blows up and have negative temperatures
    # if the max case is hit we divide by zero which yields inf so clearly flagged in any case
    @. u_or_v10 = u_or_v_bottom .* log(10 / z₀) ./ log.(max.(z_bottom, z₀) / z₀)

    # interpolate onto the output grid on the model's architecture
    return interpolate_output!(output, field, u_or_v10)
end

# add the height of the bottom layer k above the surface, T_k * R/g * log(pₛ/p_k), to z_bottom
@kernel inbounds = true function bottom_layer_height_kernel!(z_bottom, temperature, surface_pressure, coordinate, k, R_dry_over_g)
    ij = @index(Global, Linear)
    pₛ = surface_pressure[ij]
    z_bottom[ij] += temperature[ij, k] * R_dry_over_g * log(pₛ / pressure(k, pₛ, coordinate))
end

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceTemperature{F} <: AbstractOutputVariable
    name::String = "tsurf"
    unit::String = "degC"
    long_name::String = "Surface air temperature"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 12
    transform::F = (x) -> x - oftype(x, 273.15)     # [K] to [˚C]
end

# not the actual surface temperature but the core variable to read in
path(::SurfaceTemperature, simulation) = simulation.variables.grid.temperature

function interpolate_output!(
        field::AbstractField,
        output::AbstractOutput,
        variable::SurfaceTemperature,
        simulation::AbstractSimulation,
    )
    # reuse scratch array to avoid allocations
    Ts = simulation.variables.scratch.grid.a_2D

    # Retrieve T_bottom
    T_grid = get_prognostic_step(path(variable, simulation), simulation.model.time_stepping, output)
    nlayers = size(T_grid, 2)

    # Other parameters
    κ = simulation.model.atmosphere.κ
    coord = simulation.model.geometry.vertical_coordinates
    pₛ = simulation.variables.parameterizations.surface_pressure

    # Compute Ts assuming dry adiabatic profile: T_surf = T_bottom * (pₛ / p_bottom)^κ
    # in a kernel so that the layer index stays an integer (Reactant traces numbers in broadcasts)
    launch!(
        architecture(Ts), LinearWorkOrder, size(Ts), surface_air_temperature_kernel!,
        Ts, T_grid, pₛ, coord, nlayers, κ
    )
    (; transform) = variable
    @. Ts = transform(Ts)   # Convert to °C

    # interpolate onto the output grid on the model's architecture
    return interpolate_output!(output, field, Ts)
end

# dry-adiabatic descent from layer k to the surface, T_k * (pₛ/p_k)^κ, at every grid point ij
@kernel inbounds = true function surface_air_temperature_kernel!(Ts, temperature, surface_pressure, coordinate, k, κ)
    ij = @index(Global, Linear)
    pₛ = surface_pressure[ij]
    Ts[ij] = temperature[ij, k] * (pₛ / pressure(k, pₛ, coordinate))^κ
end


"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MomentumBoundaryLayerDrag <: AbstractOutputVariable
    name::String = "bld_m"
    unit::String = "1"
    long_name::String = "Momentum flux related boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::MomentumBoundaryLayerDrag, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag_momentum

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct HeatBoundaryLayerDrag <: AbstractOutputVariable
    name::String = "bld_h"
    unit::String = "1"
    long_name::String = "Heat flux related boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::HeatBoundaryLayerDrag, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag_heat

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct HumidityBoundaryLayerDrag <: AbstractOutputVariable
    name::String = "bld_q"
    unit::String = "1"
    long_name::String = "Moisture flux related boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::HumidityBoundaryLayerDrag, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag_humidity

"""Defines netCDF output of the boundary layer height above the surface, diagnosed by the
vertical diffusion (`BulkRichardsonDiffusion`) where the bulk Richardson number relative to the
lowermost layer reaches its critical value. See [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct BoundaryLayerHeight <: AbstractOutputVariable
    name::String = "blh"
    unit::String = "m"
    long_name::String = "Boundary layer height above the surface"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::BoundaryLayerHeight, simulation) =
    simulation.variables.parameterizations.boundary_layer_height

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct MomentumSurfaceRoughness <: AbstractOutputVariable
    name::String = "sr_m"
    unit::String = "m"
    long_name::String = "Momentum surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::MomentumSurfaceRoughness, simulation) =
    simulation.variables.parameterizations.momentum_roughness

@kwdef mutable struct HeatSurfaceRoughness <: AbstractOutputVariable
    name::String = "sr_h"
    unit::String = "m"
    long_name::String = "Heat surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::HeatSurfaceRoughness, simulation) =
    simulation.variables.parameterizations.heat_roughness

# collect in groups for convenience
SurfaceRoughness() = (
    MomentumSurfaceRoughness(),
    HeatSurfaceRoughness(),
)

BoundaryLayerDrag() = (
    MomentumBoundaryLayerDrag(),
    HeatBoundaryLayerDrag(),
    HumidityBoundaryLayerDrag(),
)

# collect all in one for convenience
BoundaryLayer() = (
    ZonalVelocity10m(),
    MeridionalVelocity10m(),
    SurfaceTemperature(),
    BoundaryLayerHeight(),
    BoundaryLayerDrag()...,
    SurfaceRoughness()...,
)
