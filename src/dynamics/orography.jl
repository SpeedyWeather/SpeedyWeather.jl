"""Orography with zero height in `orography` and zero surface geopotential `geopot_surf`.
$(TYPEDFIELDS)"""
struct NoOrography{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF,Grid}
    "height [m] on grid-point space."
    orography::Grid
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::LowerTriangularMatrix{Complex{NF}}
end

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function NoOrography(spectral_grid::SpectralGrid)
    (;NF, Grid, nlat_half, trunc) = spectral_grid
    orography   = zeros(Grid{NF},nlat_half)
    geopot_surf = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    return NoOrography{NF,Grid{NF}}(orography,geopot_surf)
end

function Base.show(io::IO,orog::AbstractOrography)
    println(io,"$(typeof(orog)) <: AbstractOrography")
    keys = propertynames(orog)
    print_fields(io,orog,keys)
end

# no further initialization needed
initialize!(::NoOrography,::AbstractPlanet,::SpectralTransform,::Geometry) = nothing


"""Zonal ridge orography after Jablonowski and Williamson, 2006.
$(TYPEDFIELDS)"""
Base.@kwdef struct ZonalRidge{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF,Grid}
    
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

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function ZonalRidge(spectral_grid::SpectralGrid;kwargs...)
    (;NF, Grid, nlat_half, trunc) = spectral_grid
    orography   = zeros(Grid{NF},nlat_half)
    geopot_surf = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    return ZonalRidge{NF,Grid{NF}}(;orography,geopot_surf,kwargs...)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`,`geopot_surf` in `orog` following 
Jablonowski and Williamson, 2006.
"""
function initialize!(   orog::ZonalRidge,
                        P::AbstractPlanet,
                        S::SpectralTransform,
                        G::Geometry)
    
    (;gravity, rotation) = P
    (;radius) = G
    φ = G.latds                         # latitude for each grid point [˚N]

    (;orography, geopot_surf, η₀, u₀) = orog

    ηᵥ = (1-η₀)*π/2                     # ηᵥ-coordinate of the surface [1]
    A = u₀*cos(ηᵥ)^(3/2)                # amplitude [m/s]
    RΩ = radius*rotation                # [m/s]
    g⁻¹ = inv(gravity)                  # inverse gravity [s²/m]

    for ij in eachindex(φ,orography)
        sinφ = sind(φ[ij])
        cosφ = cosd(φ[ij])

        # Jablonowski & Williamson, 2006, eq. (7)
        orography[ij] = g⁻¹*A*(A*(-2*sinφ^6*(cosφ^2 + 1/3) + 10/63) + (8/5*cosφ^3*(sinφ^2 + 2/3) - π/4)*RΩ)
    end

    spectral!(geopot_surf,orography,S)      # to grid-point space
    geopot_surf .*= gravity                 # turn orography into surface geopotential
    spectral_truncation!(geopot_surf)       # set the lmax+1 harmonics to zero
end


"""Earth's orography read from file, with smoothing.
$(TYPEDFIELDS)"""
Base.@kwdef struct EarthOrography{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF,Grid}

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

    "resolution of orography in spectral trunc"
    smoothing_truncation::Int = 85

    # FIELDS (to be initialized in initialize!)
    "height [m] on grid-point space."
    orography::Grid
    
    "surface geopotential, height*gravity [m²/s²]"
    geopot_surf::LowerTriangularMatrix{Complex{NF}} 
end

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function EarthOrography(spectral_grid::SpectralGrid;kwargs...)
    (;NF, Grid, nlat_half, trunc) = spectral_grid
    orography   = zeros(Grid{NF},nlat_half)
    geopot_surf = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    return EarthOrography{NF,Grid{NF}}(;orography,geopot_surf,kwargs...)
end

"""
$(TYPEDSIGNATURES)
Initialize the arrays `orography`,`geopot_surf` in `orog` by reading the
orography field from file.
"""
function initialize!(   orog::EarthOrography,
                        P::AbstractPlanet,
                        S::SpectralTransform,
                        G::Geometry)

    (;orography, geopot_surf) = orog
    (;gravity) = P

    # LOAD NETCDF FILE
    if orog.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__,"../../input_data",orog.file)
    else
        path = joinpath(orog.path,orog.file)
    end
    ncfile = NCDataset(path)
    orography_highres = ncfile["orog"][:,:]        # height [m]

    # Interpolate/coarsen to desired resolution
    # TODO also read lat,lon from file and flip array in case it's not as expected
    recompute_legendre = true   # don't allocate large arrays as spectral transform is not reused
    orography_spec = spectral(orography_highres;Grid=orog.file_Grid,recompute_legendre)
    
    copyto!(geopot_surf,orography_spec)     # truncates to the size of geopot_surf, no *gravity yet
    if orog.smoothing                       # smooth orography in spectral space?
        SpeedyTransforms.spectral_smoothing!(geopot_surf,orog.smoothing_strength,
                                                            power=orog.smoothing_power,
                                                            truncation=orog.smoothing_truncation)
    end

    gridded!(orography,geopot_surf,S)       # to grid-point space
    geopot_surf .*= gravity                 # turn orography into surface geopotential
    spectral_truncation!(geopot_surf)       # set the lmax+1 harmonics to zero
end

"""Land-sea mask, fractional, read from file.
$(TYPEDFIELDS)"""
Base.@kwdef struct LandSeaMask{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractLandSeaMask{NF,Grid}

    # OPTIONS
    "path to the folder containing the land-sea mask file, pkg path default"
    path::String = "SpeedyWeather.jl/input_data"

    "filename of land sea mask"
    file::String = "land-sea_mask.nc"

    "Grid the land-sea mask file comes on"
    file_Grid::Type{<:AbstractGrid} = FullClenshawGrid

    # FIELDS (to be initialized in initialize!)
    "Land-sea mask [1] on grid-point space. Land=1, sea=0, land-area fraction in between."
    land_sea_mask::Grid
end

"""
$(TYPEDSIGNATURES)
Generator function pulling the resolution information from `spectral_grid`."""
function LandSeaMask(spectral_grid::SpectralGrid;kwargs...)
    (;NF, Grid, nlat_half) = spectral_grid
    land_sea_mask   = zeros(Grid{NF},nlat_half)
    return LandSeaMask{NF,Grid{NF}}(;land_sea_mask,kwargs...)
end

function initialize!(land_sea_mask::LandSeaMask)

    (;file_Grid) = land_sea_mask

    # LOAD NETCDF FILE
    if land_sea_mask.path == "SpeedyWeather.jl/input_data"
        path = joinpath(@__DIR__,"../../input_data",land_sea_mask.file)
    else
        path = joinpath(land_sea_mask.path,land_sea_mask.file)
    end
    ncfile = NCDataset(path)
    
    # high resolution land-sea mask
    lsm_highres = file_Grid(ncfile["lsm"][:,:])

    # average onto grid cells of the model
    RingGrids.grid_cell_average!(land_sea_mask.land_sea_mask,lsm_highres)
end