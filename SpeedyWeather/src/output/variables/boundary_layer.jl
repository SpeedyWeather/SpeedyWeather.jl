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

    # Compute z_bottom as z_surf + T_bottom * Δp_geopot / g, start with z_surf
    z_bottom .= max.(simulation.model.orography.orography, 0)   # [m] set negative values to zero
    T = get_prognostic_step(simulation.variables.grid.temperature, simulation.model.time_stepping, output) 
    T_bottom = field_view(T, :, nlayers)
    # TODO: Δp_geopot_full[end] = R*log(σ_half/σ_full) at the bottom layer is a precomputed
    # sigma-only constant. Generalising to hybrid coordinates requires computing R*log(p_half/p_full)
    # at each grid point using the surface pressure field.
    Δp_geopot = simulation.model.geopotential.Δp_geopot_full[end]

    # accumulate the second term in
    @. z_bottom += T_bottom * Δp_geopot / simulation.model.planet.gravity

    # Compute u10, TODO should this be the same z₀ as in vertical diffusion or surface fluxes?
    z₀ = simulation.model.vertical_diffusion.roughness_length

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
@kwdef mutable struct BoundaryLayerDragOutput <: AbstractOutputVariable
    name::String = "bld"
    unit::String = "1"
    long_name::String = "Boundary layer drag coefficient"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::BoundaryLayerDragOutput, simulation) =
    simulation.variables.parameterizations.boundary_layer_drag

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct LandSurfaceRoughnessOutput <: AbstractOutputVariable
    name::String = "lsr"
    unit::String = "m"
    long_name::String = "Land surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::LandSurfaceRoughnessOutput, simulation) =
    simulation.variables.parameterizations.land.surface_roughness

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct OceanSurfaceRoughnessOutput <: AbstractOutputVariable
    name::String = "osr"
    unit::String = "m"
    long_name::String = "Ocean surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::OceanSurfaceRoughnessOutput, simulation) =
    simulation.variables.parameterizations.ocean.surface_roughness

"""Defines netCDF output for a specific variables, see [`VorticityOutput`](@ref) for details.
Fields are: $(TYPEDFIELDS)"""
@kwdef mutable struct SurfaceRoughnessOutput <: AbstractOutputVariable
    name::String = "sr"
    unit::String = "m"
    long_name::String = "Surface roughness length"
    dims_xyzt::NTuple{4, Bool} = (true, true, false, true)
    missing_value::Float64 = NaN
    compression_level::Int = 3
    shuffle::Bool = true
    keepbits::Int = 7
end

path(::SurfaceRoughnessOutput, simulation) =
    simulation.variables.parameterizations.surface_roughness

# collect all in one for convenience
BoundaryLayerOutput() = (
    ZonalVelocity10mOutput(),
    MeridionalVelocity10mOutput(),
    SurfaceTemperatureOutput(),
    BoundaryLayerDragOutput(),
    SurfaceRoughnessOutput(),
)
