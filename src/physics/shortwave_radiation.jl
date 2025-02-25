abstract type AbstractShortwave <: AbstractRadiation end

function get_nbands(R::AbstractRadiation)
    hasfield(typeof(R), :nbands) && return R.nbands
    return 0
end

# function barrier for all AbstractShortwave
function shortwave_radiation!(column::ColumnVariables, model::PrimitiveEquation)
    shortwave_radiation!(column, model.shortwave_radiation, model)
end

## NO SHORTWAVE RADIATION
export NoShortwave
struct NoShortwave <: AbstractShortwave end
NoShortwave(SG::SpectralGrid) = NoShortwave()
initialize!(::NoShortwave, ::PrimitiveEquation) = nothing
shortwave_radiation!(::ColumnVariables, ::NoShortwave, ::PrimitiveEquation) = nothing

## SHORTWAVE RADIATION FOR A FULLY TRANSPARENT ATMOSPHERE
export TransparentShortwave
struct TransparentShortwave <: AbstractShortwave end
TransparentShortwave(SG::SpectralGrid) = TransparentShortwave()
initialize!(::TransparentShortwave, ::PrimitiveEquation) = nothing

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    model::PrimitiveEquation,
)
    shortwave_radiation!(column, scheme, model.planet)
end

function shortwave_radiation!(
    column::ColumnVariables,
    scheme::TransparentShortwave,
    planet::AbstractPlanet,
)
    (; cos_zenith, land_fraction, albedo_ocean, albedo_land) = column
    (; solar_constant) = planet

    # ocean
    column.surface_shortwave_down_ocean = (1 - albedo_ocean) * solar_constant * cos_zenith
    column.surface_shortwave_up_ocean = albedo_ocean * solar_constant * cos_zenith

    # land
    column.surface_shortwave_down_land = (1 - albedo_land) * solar_constant * cos_zenith
    column.surface_shortwave_up_land = albedo_land * solar_constant * cos_zenith

    # land-sea mask-weighted OSR
    column.outgoing_shortwave_radiation = (1 - land_fraction)*column.surface_shortwave_up_ocean +
                                            land_fraction*column.surface_shortwave_up_land
    return nothing
end

# NBandRadiation is defined in longwave_radiation.jl
function shortwave_radiation!(
    column::ColumnVariables,
    scheme::NBandRadiation,
    model::PrimitiveEquation,
)

    (; nlayers, cos_zenith, albedo) = column
    nbands = column.nbands_shortwave                # number of spectral bands
    (; solar_constant) = model.planet

    @inbounds for band in 1:nbands                  # loop over spectral bands
        dτ = view(column.optical_depth_shortwave, :, band)   # differential optical depth per layer of that band

        # DOWNWARD flux D
        D = solar_constant * cos_zenith             # top boundary condition of longwave flux
        column.flux_temp_downward[1] += D           # accumulate the top downward flux

        for k in 1:nlayers
            D -= dτ[k]*D                            # flux through layer k with optical depth dτ, radiative transfer
            column.flux_temp_downward[k+1] += D
        end

        # UPWARD flux U
        U = D*albedo                                # boundary condition at surface, reflection from albedo
        column.flux_temp_upward[nlayers+1] += U     # accumulate fluxes

        for k in nlayers:-1:1                       # integrate from surface up
            # Radiative transfer, e.g. Frierson et al. 2006, equation 6
            U -= dτ[k]*U                            # negative because we integrate from surface up in -τ direction
            column.flux_temp_upward[k] += U         # accumulate that flux
        end

        # store outgoing shortwave radiation (OSR) for diagnostics, accumulate over bands (reset when column is reset)
        column.outgoing_shortwave_radiation += U
    end
end