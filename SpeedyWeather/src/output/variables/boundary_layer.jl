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
function output!(
        output::AbstractOutput,
        variable::Union{ZonalVelocity10m, MeridionalVelocity10m},
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
    z₀ = simulation.variables.parameterizations.surface_roughness

    # include a max here as this will throw an error if the model blows up and have negative temperatures
    # if the max case is hit we divide by zero which yields inf so clearly flagged in any case
    @. u_or_v10 = u_or_v_bottom .* log(10 / z₀) ./ log.(max.(z_bottom, z₀) / z₀)

    # interpolate onto the output grid on the model's architecture, then to host, round and write
    interpolate_output!(output, output.field2D, u_or_v10)
    write_output!(output, variable, output.field2D, output.host2D)
    return nothing
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

function output!(
        output::AbstractOutput,
        variable::SurfaceTemperature,
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

    # interpolate onto the output grid on the model's architecture, then to host, round and write
    interpolate_output!(output, output.field2D, Ts)
    write_output!(output, variable, output.field2D, output.host2D)
    return nothing
end


"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct BoundaryLayerDrag <: AbstractOutputVariable
    name::String = "bld"
    unit::String = "1"
    long_name::String = "Boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::BoundaryLayerDrag, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LandSurfaceRoughness <: AbstractOutputVariable
    name::String = "lsr"
    unit::String = "m"
    long_name::String = "Land surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::LandSurfaceRoughness, simulation) =
    simulation.variables.parameterizations.land.surface_roughness

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct OceanSurfaceRoughness <: AbstractOutputVariable
    name::String = "osr"
    unit::String = "m"
    long_name::String = "Ocean surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::OceanSurfaceRoughness, simulation) =
    simulation.variables.parameterizations.ocean.surface_roughness

"""Defines netCDF output for a specific variables, see [`Vorticity`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceRoughness <: AbstractOutputVariable
    name::String = "sr"
    unit::String = "m"
    long_name::String = "Surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceRoughness, simulation) =
    simulation.variables.parameterizations.surface_roughness

# collect all in one for convenience
BoundaryLayer() = (
    ZonalVelocity10m(),
    MeridionalVelocity10m(),
    SurfaceTemperature(),
    BoundaryLayerDrag(),
    SurfaceRoughness(),
)
