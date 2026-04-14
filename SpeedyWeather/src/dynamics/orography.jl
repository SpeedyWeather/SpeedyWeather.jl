abstract type AbstractOrography <: AbstractModelComponent end

# adapt to GPU only the fields themselves
Adapt.adapt_structure(to, orog::AbstractOrography) = (orography = adapt_structure(to, orog.orography), surface_geopotential = adapt_structure(to, orog.surface_geopotential))

# general constructor with empty arrays (will be initialized in initialize!)
function (O::Type{<:AbstractOrography})(spectral_grid::SpectralGrid; kwargs...)
    (; GridVariable2D, SpectralVariable2D, grid, spectrum) = spectral_grid
    orography = zeros(GridVariable2D, grid)
    surface_geopotential = zeros(SpectralVariable2D, spectrum)
    return O(; orography, surface_geopotential, kwargs...)
end

export NoOrography

"""Orography with zero height in `orography` and zero surface geopotential `surface_geopotential`.
$(TYPEDFIELDS)"""
@kwdef struct NoOrography{G, S} <: AbstractOrography
    "[OPTION] height [m] on grid-point space."
    orography::G

    "[OPTION] surface geopotential, height*gravity [m²/s²]"
    surface_geopotential::S
end

# set arrays to zero when initialized
function initialize!(O::NoOrography, ::AbstractModel)
    O.orography .= 0
    O.surface_geopotential .= 0
    return nothing
end

export ManualOrography

"""Orography with zero height in `orography` and zero surface geopotential `surface_geopotential`.
$(TYPEDFIELDS)"""
@kwdef struct ManualOrography{G, S} <: AbstractOrography
    "[OPTION] height [m] on grid-point space."
    orography::G

    "[OPTION] surface geopotential, height*gravity [m²/s²]"
    surface_geopotential::S
end

# deliberate don't touch the arrays in initialize, as they will be set by the user with set! before the model is run
initialize!(orog::ManualOrography, model::AbstractModel) = nothing

# set orography with grid, scalar, function
function set!(
        orography::AbstractOrography,
        v,                                      # new orography, function, scalar, grid, ...
        geometry::Geometry,
        spectral_transform::SpectralTransform;
        gravity = DEFAULT_GRAVITY,
        add::Bool = false,
    )
    set!(orography.orography, v, geometry, spectral_transform; add)

    # now synchronize with geopotential in spectral space (used in primitive models)
    transform!(orography.surface_geopotential, orography.orography, spectral_transform)
    orography.surface_geopotential .*= gravity
    SpeedyTransforms.spectral_truncation!(orography.surface_geopotential)     # set the lmax+1 harmonics to zero
    return nothing
end

export ZonalRidge

"""Zonal ridge orography after Jablonowski and Williamson, 2006.
$(TYPEDFIELDS)"""
@kwdef struct ZonalRidge{NF, GridVariable2D, SpectralVariable2D} <: AbstractOrography

    "[OPTION] conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::NF = 0.252

    "[OPTION] max amplitude of zonal wind [m/s] that scales orography height"
    u₀::NF = 35

    # FIELDS (to be initialized in initialize!)
    "height [m] on grid-point space."
    orography::GridVariable2D

    "surface geopotential, height*gravity [m²/s²]"
    surface_geopotential::SpectralVariable2D
end

function ZonalRidge(spectral_grid::SpectralGrid; kwargs...)
    (; architecture, NF, GridVariable2D, SpectralVariable2D, grid, spectrum) = spectral_grid
    orography = zeros(GridVariable2D, grid)
    surface_geopotential = zeros(SpectralVariable2D, spectrum)
    return ZonalRidge{NF, GridVariable2D, SpectralVariable2D}(; orography, surface_geopotential, kwargs...)
end

# function barrier
function initialize!(
        orog::ZonalRidge,
        model::AbstractModel
    )
    return initialize!(orog, model.planet, model.spectral_transform, model.geometry)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`, `surface_geopotential` in `orog` following 
Jablonowski and Williamson, 2006."""
function initialize!(
        orog::ZonalRidge,
        P::AbstractPlanet,
        S::SpectralTransform,
        G::Geometry
    )

    (; radius, gravity, rotation) = P
    φ = G.latds                         # latitude for each grid point [˚N]

    (; orography, surface_geopotential, η₀, u₀) = orog

    ηᵥ = (1 - η₀) * π / 2                     # ηᵥ-coordinate of the surface [1]
    A = u₀ * cos(ηᵥ)^(3 / 2)                # amplitude [m/s]
    RΩ = radius * rotation                # [m/s]
    g⁻¹ = inv(gravity)                  # inverse gravity [s²/m]

    # TODO use set! to write this
    for ij in eachindex(φ, orography)
        sinφ = sind(φ[ij])
        cosφ = cosd(φ[ij])

        # Jablonowski & Williamson, 2006, eq. (7)
        orography[ij] = g⁻¹ * A * (A * (-2 * sinφ^6 * (cosφ^2 + 1 / 3) + 10 / 63) + (8 / 5 * cosφ^3 * (sinφ^2 + 2 / 3) - π / 4) * RΩ)
    end

    transform!(surface_geopotential, orography, S)   # to grid-point space
    surface_geopotential .*= gravity                 # turn orography into surface geopotential
    SpeedyTransforms.spectral_truncation!(surface_geopotential)       # set the lmax+1 harmonics to zero
    return nothing
end

export EarthOrography

"""Earth's orography read from file, with smoothing.
$(TYPEDFIELDS)"""
@kwdef struct EarthOrography{NF, GridVariable2D, SpectralVariable2D} <: AbstractOrography
    "[OPTION] filename of orography"
    file::String = "orography.nc"

    "[OPTION] path to the folder containing the orography"
    path::String = joinpath("data", "boundary_conditions", file)

    "[OPTION] flag to check for orography in SpeedyWeatherAssets or locally"
    from_assets::Bool = true

    "[OPTION] SpeedyWeatherAssets version number"
    version::VersionNumber = DEFAULT_ASSETS_VERSION

    "[OPTION] NCDataset variable name"
    varname::String = "orog"

    "[OPTION] Grid the orography file comes on"
    FieldType::Type{<:AbstractField} = FullGaussianField

    "[OPTION] scale orography by a factor"
    scale::NF = 1.0

    "[OPTION] smooth the orography field?"
    smoothing::Bool = true

    "[OPTION] ower of Laplacian for smoothing"
    smoothing_power::NF = 1.0

    "[OPTION] highest degree l is multiplied by"
    smoothing_strength::NF = 0.1

    "[OPTION] fraction of highest wavenumbers to smooth"
    smoothing_fraction::NF = 0.05

    # FIELDS (to be initialized in initialize!)
    "[DERIVED] height [m] on grid-point space."
    orography::GridVariable2D

    "[DERIVED] surface geopotential, height*gravity [m²/s²]"
    surface_geopotential::SpectralVariable2D
end

function EarthOrography(spectral_grid::SpectralGrid; kwargs...)
    (; NF, GridVariable2D, SpectralVariable2D, grid, spectrum) = spectral_grid
    orography = zeros(GridVariable2D, grid)
    surface_geopotential = zeros(SpectralVariable2D, spectrum)
    return EarthOrography{NF, GridVariable2D, SpectralVariable2D}(; orography, surface_geopotential, kwargs...)
end

# function barrier
function initialize!(
        orog::EarthOrography,
        model::AbstractModel
    )
    return initialize!(orog, model.planet, model.spectral_transform)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`, `surface_geopotential` in `orog` by reading the
orography field from file.
"""
function initialize!(
        orog::EarthOrography,
        P::AbstractPlanet,
        S::AbstractSpectralTransform
    )

    (; orography, surface_geopotential, scale) = orog
    (; gravity) = P

    # load orography from file
    field = get_asset(
        orog.path;
        from_assets = orog.from_assets,
        name = orog.varname,
        ArrayType = orog.FieldType,
        FileFormat = NCDataset,
        version = orog.version
    )

    orography_highres = on_architecture(S.architecture, field)

    # Interpolate/coarsen to desired resolution
    interpolate!(orography, orography_highres)
    orography .*= scale                     # scale orography (default 1)
    transform!(surface_geopotential, orography, S)   # no *gravity yet

    if orog.smoothing                       # smooth orography in spectral space?
        # get trunc=lmax from size of surface_geopotential
        trunc = (size(surface_geopotential, 1, as = Matrix) - 2)
        # degree of harmonics to be truncated
        truncation = round(Int, trunc * (1 - orog.smoothing_fraction))
        c = orog.smoothing_strength
        power = orog.smoothing_power
        SpeedyTransforms.spectral_smoothing!(surface_geopotential, c; power, truncation)
    end

    transform!(orography, surface_geopotential, S)                  # to grid-point space
    surface_geopotential .*= gravity                                # turn orography into surface geopotential
    SpeedyTransforms.spectral_truncation!(surface_geopotential)     # set the lmax+1 harmonics to zero
    return nothing
end
