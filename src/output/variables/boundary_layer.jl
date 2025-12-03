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

path(::ZonalVelocity10mOutput, simulation) = simulation.diagnostic_variables.grid.u_grid

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

path(::MeridionalVelocity10mOutput, simulation) = simulation.diagnostic_variables.grid.v_grid

"""$(TYPEDSIGNATURES)
10m wind is defined using a logarithmic profile from the lowermost model layer.
```math
u_{10}=u_{bottom} * ln(10/z₀) / ln(z_{bottom}/z₀)
``` (with z in m), 
with 
```math
z_{bottom} = z_{surf} + T_{bottom} * Δp_geopot_{bottom} / g
```
"""
function output!(
    output::NetCDFOutput,
    variable::Union{ZonalVelocity10mOutput, MeridionalVelocity10mOutput},
    simulation::AbstractSimulation,
)
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    # Retrieve u_bottom
    u_or_v_grid = path(variable, simulation)
    nlayers = size(u_or_v_grid, 2)
    u_or_v_bottom = field_view(u_or_v_grid, :, nlayers)
    
    z_bottom = simulation.diagnostic_variables.dynamics.a_2D_grid
    u_or_v10 = simulation.diagnostic_variables.dynamics.b_2D_grid
    
    # Compute z_bottom as z_surf + T_bottom * Δp_geopot / g, start with z_surf
    z_bottom .= max.(simulation.model.orography.orography, 0)   # [m] set negative values to zero
    T_bottom = field_view(simulation.diagnostic_variables.grid.temp_grid, :, nlayers)
    Δp_geopot = simulation.model.geopotential.Δp_geopot_full[end]
    
    # accumulate the second term in
    @. z_bottom += T_bottom * Δp_geopot / simulation.model.planet.gravity

    # Compute u10, TODO should this be the same z₀ as in vertical diffusion or surface fluxes?
    z₀ = simulation.model.vertical_diffusion.roughness_length
    @. u_or_v10 = u_or_v_bottom .* log(10/z₀) ./ log.(z_bottom/z₀)

    # interpolate 2D/3D variables
    u_or_v10_output = output.field2D
    u_or_v10_grid = on_architecture(CPU(), u_or_v10)
    RingGrids.interpolate!(u_or_v10_output, u_or_v10_grid, output.interpolator)

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(u_or_v10_output, variable.keepbits)
    end

    i = output.output_counter               # output time step i to write
    indices = get_indices(i, variable)      # returns (:, :, i) for example, depending on dims
    output.netcdf_file[variable.name][indices...] = u_or_v10_output     # actually write to file
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
path(::SurfaceTemperatureOutput, simulation) = simulation.diagnostic_variables.grid.temp_grid

function output!(
    output::NetCDFOutput,
    variable::SurfaceTemperatureOutput,
    simulation::AbstractSimulation,
)
    # escape immediately after first call if variable doesn't have a time dimension
    ~hastime(variable) && output.output_counter > 1 && return nothing

    # reuse scratch array to avoid allocations
    Ts = simulation.diagnostic_variables.dynamics.a_2D_grid

    # Retrieve T_bottom
    T_grid = path(variable, simulation)
    nlayers = size(T_grid, 2)
    T_bottom = field_view(T_grid, :, nlayers)

    # Other parameters
    κ = simulation.model.atmosphere.κ
    σ_bottom = simulation.model.geometry.σ_levels_full[end]

    # Compute Ts
    (; transform) = variable
    @. Ts = transform(T_bottom * σ_bottom ^ (-κ))   # Convert to °C

    # interpolate 2D/3D variables
    Ts_output = output.field2D
    Ts_grid = on_architecture(CPU(), Ts)
    RingGrids.interpolate!(Ts_output, Ts_grid, output.interpolator)

    if hasproperty(variable, :keepbits)     # round mantissabits for compression
        round!(Ts_output, variable.keepbits)
    end

    i = output.output_counter               # output time step i to write
    indices = get_indices(i, variable)      # returns (:, :, i) for example, depending on dims
    output.netcdf_file[variable.name][indices...] = Ts_output     # actually write to file
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
    simulation.diagnostic_variables.physics.boundary_layer_drag

# collect all in one for convenience
BoundaryLayerOutput() = (
    ZonalVelocity10mOutput(),
    MeridionalVelocity10mOutput(),
    SurfaceTemperatureOutput(),
    BoundaryLayerDragOutput(),
)