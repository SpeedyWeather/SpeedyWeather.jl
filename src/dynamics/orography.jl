abstract type AbstractOrography <: AbstractModelComponent end
export NoOrography

"""Orography with zero height in `orography` and zero surface geopotential `geopot_surf`.
$(TYPEDFIELDS)"""
@kwdef struct NoOrography{NF, GridVariable2D, SpectralVariable2D} <: AbstractOrography
    "height [m] on grid-point space."
    orography::GridVariable2D
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::SpectralVariable2D
end

# constructor
function NoOrography(spectral_grid::SpectralGrid)
    (; NF, GridVariable2D, SpectralVariable2D, nlat_half, trunc) = spectral_grid
    orography   = zeros(GridVariable2D, nlat_half)
    geopot_surf = zeros(SpectralVariable2D, trunc+2, trunc+1)
    return NoOrography{NF, GridVariable2D, SpectralVariable2D}(; orography, geopot_surf)
end

# no further initialization needed
initialize!(::NoOrography, ::AbstractModel) = nothing

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
    transform!(orography.geopot_surf, orography.orography, spectral_transform)
    orography.geopot_surf .*= gravity
    spectral_truncation!(orography.geopot_surf)     # set the lmax+1 harmonics to zero
    return nothing
end

export ZonalRidge

"""Zonal ridge orography after Jablonowski and Williamson, 2006.
$(TYPEDFIELDS)"""
@kwdef struct ZonalRidge{NF, GridVariable2D, SpectralVariable2D} <: AbstractOrography
    
    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::NF = 0.252

    "max amplitude of zonal wind [m/s] that scales orography height"
    u₀::NF = 35

    # FIELDS (to be initialized in initialize!)
    "height [m] on grid-point space."
    orography::GridVariable2D
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::SpectralVariable2D
end

# constructor
function ZonalRidge(spectral_grid::SpectralGrid; kwargs...)
    (; NF, GridVariable2D, SpectralVariable2D, nlat_half, trunc) = spectral_grid
    orography   = zeros(GridVariable2D, nlat_half)
    geopot_surf = zeros(SpectralVariable2D, trunc+2, trunc+1)
    return ZonalRidge{NF, GridVariable2D, SpectralVariable2D}(;
        orography, geopot_surf, kwargs...)
end

# function barrier
function initialize!(   orog::ZonalRidge,
                        model::AbstractModel)
    initialize!(orog, model.planet, model.spectral_transform, model.geometry)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`, `geopot_surf` in `orog` following 
Jablonowski and Williamson, 2006."""
function initialize!(   orog::ZonalRidge,
                        P::AbstractPlanet,
                        S::SpectralTransform,
                        G::Geometry)
    
    (; gravity, rotation) = P
    (; radius) = G
    φ = G.latds                         # latitude for each grid point [˚N]

    (; orography, geopot_surf, η₀, u₀) = orog

    ηᵥ = (1-η₀)*π/2                     # ηᵥ-coordinate of the surface [1]
    A = u₀*cos(ηᵥ)^(3/2)                # amplitude [m/s]
    RΩ = radius*rotation                # [m/s]
    g⁻¹ = inv(gravity)                  # inverse gravity [s²/m]

    # TODO use set! to write this
    for ij in eachindex(φ, orography)
        sinφ = sind(φ[ij])
        cosφ = cosd(φ[ij])

        # Jablonowski & Williamson, 2006, eq. (7)
        orography[ij] = g⁻¹*A*(A*(-2*sinφ^6*(cosφ^2 + 1/3) + 10/63) + (8/5*cosφ^3*(sinφ^2 + 2/3) - π/4)*RΩ)
    end

    transform!(geopot_surf, orography, S)   # to grid-point space
    geopot_surf .*= gravity                 # turn orography into surface geopotential
    spectral_truncation!(geopot_surf)       # set the lmax+1 harmonics to zero
    return nothing
end

export EarthOrography

"""Earth's orography read from file, with smoothing.
$(TYPEDFIELDS)"""
@kwdef struct EarthOrography{NF, GridVariable2D, SpectralVariable2D} <: AbstractOrography

    # OPTIONS
    "path to the folder containing the orography file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of orography"
    file::String = "orography.nc"

    "Grid the orography file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "scale orography by a factor"
    scale::NF = 1.0

    "smooth the orography field?"
    smoothing::Bool = true

    "power of Laplacian for smoothing"
    smoothing_power::NF = 1.0

    "highest degree l is multiplied by"
    smoothing_strength::NF = 0.1

    "fraction of highest wavenumbers to smooth"
    smoothing_fraction::NF = 0.05

    # FIELDS (to be initialized in initialize!)
    "height [m] on grid-point space."
    orography::GridVariable2D
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::SpectralVariable2D
end

# constructor
function EarthOrography(spectral_grid::SpectralGrid; kwargs...)
    (; architecture, NF, GridVariable2D, SpectralVariable2D, grid, spectrum) = spectral_grid
    orography   = on_architecture(architecture, zeros(GridVariable2D, grid))
    geopot_surf = on_architecture(architecture, zeros(SpectralVariable2D, spectrum))
    return EarthOrography{NF, GridVariable2D, SpectralVariable2D}(;
        orography, geopot_surf, kwargs...)
end

# function barrier
function initialize!(   orog::EarthOrography,
                        model::AbstractModel)
    initialize!(orog, model.planet, model.spectral_transform)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`, `geopot_surf` in `orog` by reading the
orography field from file.
"""
function initialize!(   orog::EarthOrography,
                        P::AbstractPlanet,
                        S::SpectralTransform)

    (; orography, geopot_surf, scale) = orog
    (; gravity) = P

    # LOAD NETCDF FILE
    if orog.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__, "../../input_data", orog.file)
    else
        path = joinpath(orog.path, orog.file)
    end
    ncfile = NCDataset(path)

    # height [m], wrap matrix into a grid
    # TODO also read lat, lon from file and flip array in case it's not as expected
    # F = RingGrids.field_type(orog.file_Grid)  # TODO this isn't working, hardcode instead
    orography_highres = on_architecture(S.architecture, FullGaussianField(ncfile["orog"].var[:, :], input_as=Matrix))

    # Interpolate/coarsen to desired resolution
    interpolate!(orography, orography_highres)
    orography .*= scale                     # scale orography (default 1)
    transform!(geopot_surf, orography, S)   # no *gravity yet
  
    if orog.smoothing                       # smooth orography in spectral space?
        # get trunc=lmax from size of geopot_surf
        trunc = (size(geopot_surf, 1, as=Matrix) - 2)
        # degree of harmonics to be truncated
        truncation = round(Int, trunc * (1-orog.smoothing_fraction))
        c = orog.smoothing_strength
        power = orog.smoothing_power
        SpeedyTransforms.spectral_smoothing!(geopot_surf, c; power, truncation)
    end

    transform!(orography, geopot_surf, S)   # to grid-point space
    geopot_surf .*= gravity                 # turn orography into surface geopotential
    spectral_truncation!(geopot_surf)       # set the lmax+1 harmonics to zero
    return nothing    
end