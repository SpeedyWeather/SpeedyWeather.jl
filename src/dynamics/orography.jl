abstract type AbstractOrography <: AbstractModelComponent end
export NoOrography

"""Orography with zero height in `orography` and zero surface geopotential `geopot_surf`.
$(TYPEDFIELDS)"""
@kwdef struct NoOrography{NF<:AbstractFloat, Grid<:AbstractGrid{NF}} <: AbstractOrography
    "height [m] on grid-point space."
    orography::Grid
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::LowerTriangularMatrix{Complex{NF}}
end

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid` for
all Orography <: AbstractOrography."""
function (::Type{Orography})(
    spectral_grid::SpectralGrid;
    kwargs...
) where Orography <: AbstractOrography
    (; NF, Grid, nlat_half, trunc) = spectral_grid
    orography   = zeros(Grid{NF}, nlat_half)
    geopot_surf = zeros(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)
    return Orography{NF, Grid{NF}}(; orography, geopot_surf, kwargs...)
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

# convenient function barrier to use geometry, spectral_transform and gravity from model
function set!(model::AbstractModel; orography, kwargs...)
    set!(model.orography, orography, model.geometry, model.spectral_transform;
        gravity=model.planet.gravity, kwargs...)
end

export ZonalRidge

"""Zonal ridge orography after Jablonowski and Williamson, 2006.
$(TYPEDFIELDS)"""
@kwdef struct ZonalRidge{NF<:AbstractFloat, Grid<:AbstractGrid{NF}} <: AbstractOrography
    
    "conversion from σ to Jablonowski's ηᵥ-coordinates"
    η₀::Float64 = 0.252

    "max amplitude of zonal wind [m/s] that scales orography height"
    u₀::Float64 = 35

    # FIELDS (to be initialized in initialize!)
    "height [m] on grid-point space."
    orography::Grid
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::LowerTriangularMatrix{Complex{NF}} 
end

# function barrier
function initialize!(   orog::ZonalRidge,
                        model::AbstractModel)
    initialize!(orog, model.planet, model.spectral_transform, model.geometry)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`, `geopot_surf` in `orog` following 
Jablonowski and Williamson, 2006.
"""
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
@kwdef struct EarthOrography{NF<:AbstractFloat, Grid<:AbstractGrid{NF}} <: AbstractOrography

    # OPTIONS
    "path to the folder containing the orography file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of orography"
    file::String = "orography.nc"

    "Grid the orography file comes on"
    file_Grid::Type{<:AbstractGrid} = FullGaussianGrid

    "scale orography by a factor"
    scale::Float64 = 1

    "smooth the orography field?"
    smoothing::Bool = true

    "power of Laplacian for smoothing"
    smoothing_power::Float64 = 1.0

    "highest degree l is multiplied by"
    smoothing_strength::Float64 = 0.1

    "fraction of highest wavenumbers to smooth"
    smoothing_fraction::Float64 = 0.05

    # FIELDS (to be initialized in initialize!)
    "height [m] on grid-point space."
    orography::Grid
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::LowerTriangularMatrix{Complex{NF}} 
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
    orography_highres = orog.file_Grid(ncfile["orog"].var[:, :], input_as=Matrix)

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