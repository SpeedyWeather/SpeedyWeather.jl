abstract type AbstractShortwave <: AbstractRadiation end

function get_nbands(radiation::Union{AbstractRadiation, Nothing})
    hasfield(typeof(radiation), :nbands) && return radiation.nbands
    return 0
end

## SHORTWAVE RADIATION FOR A FULLY TRANSPARENT ATMOSPHERE
export TransparentShortwave
struct TransparentShortwave <: AbstractShortwave end
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()

function variables(::TransparentShortwave)
    return (
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down", units="W/m^2"),
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down over ocean", units="W/m^2", namespace=:ocean),
        DiagnosticVariable(name=:surface_shortwave_down, dims=Grid2D(), desc="Surface shortwave radiation down over land", units="W/m^2", namespace=:land),
        DiagnosticVariable(name=:surface_shortwave_up,   dims=Grid2D(), desc="Surface shortwave radiation up",   units="W/m^2"),
        DiagnosticVariable(name=:outgoing_shortwave,     dims=Grid2D(), desc="TOA Shortwave radiation up",       units="W/m^2"),
        DiagnosticVariable(name=:cos_zenith,             dims=Grid2D(), desc="Cos zenith angle",                 units="1"),
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo over ocean", units="1", namespace=:ocean),
        DiagnosticVariable(name=:albedo, dims=Grid2D(), desc="Albedo over land", units="1", namespace=:land),
    )
end

initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

# function barrier
parameterization!(ij, diagn, progn, shortwave::TransparentShortwave, model) =
    shortwave_radiation!(ij, diagn, shortwave, model.planet, model.land_sea_mask)

function shortwave_radiation!(ij, diagn, shortwave::TransparentShortwave, planet, land_sea_mask)

    (; surface_shortwave_down, surface_shortwave_up) = diagn.physics
    (; outgoing_shortwave) = diagn.physics
    ssrd_ocean = diagn.physics.ocean.surface_shortwave_down
    ssrd_land = diagn.physics.land.surface_shortwave_down

    cos_zenith = diagn.physics.cos_zenith[ij]
    land_fraction = land_sea_mask[ij]
    albedo_ocean = diagn.physics.ocean.albedo[ij]
    albedo_land = diagn.physics.land.albedo[ij]
    S₀ = planet.solar_constant

    D = S₀ * cos_zenith             # top of atmosphere downward radiation
    surface_shortwave_down[ij] = D  # transparent atmosphere so same at surface (before albedo)

    # shortwave up is after albedo reflection, separated by ocean/land
    ssrd_ocean[ij] = albedo_ocean * D
    ssrd_land[ij]  = albedo_land * D

    # land-sea mask-weighted
    albedo = (1 - land_fraction)*albedo_ocean + land_fraction*albedo_land
    surface_shortwave_up[ij] = albedo * D

    # transparent also for reflected shortwave radiation travelling up
    outgoing_shortwave[ij] = surface_shortwave_up[ij]
    return nothing
end