abstract type AbstractOrography{NF} end

# gridded field (orography) and spectral field (surface geopotential) for Earth's orography
struct EarthOrography{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF}
    orography::Grid                                 # height [m]
    geopot_surf::LowerTriangularMatrix{Complex{NF}} # surface geopotential, height*gravity [m²/s²]
end

abstract type ZonalRidge{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF} end
ZonalRidge(orography,geopot_surf) = EarthOrography(orography,geopot_surf)

abstract type NoOrography{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF} end
NoOrography(orography,geopot_surf) = EarthOrography(orography,geopot_surf)

function zero(O::Type{<:AbstractOrography},S::SpectralTransform{NF}) where NF
    @unpack Grid, nlat_half, lmax, mmax = S
    orography   = zeros(Grid{NF},nlat_half)
    geopot_surf = zeros(LowerTriangularMatrix{Complex{NF}},lmax+1,mmax)
    return O(orography,geopot_surf)
end

abstract type AbstractBoundaries{NF} end

struct Boundaries{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractBoundaries{NF}
    orography::AbstractOrography{NF}
    # landsea_mask::AbstractLandSeaMask{NF}
    # albedo::AbstractAlbedo{NF}
end

function Boundaries(P::Parameters,
                    S::SpectralTransform,
                    G::Geometry)

    orography = zero(P.orography,S)
    initialize_orography!(orography,P.orography,P,S,G)

    return Boundaries(orography)
end

function initialize_orography!( ::AbstractOrography,
                                ::Type{<:NoOrography},
                                ::Parameters{<:ModelSetup},
                                args...)
    return nothing
end

function initialize_orography!( ::AbstractOrography,
                                ::Type{<:AbstractOrography},
                                ::Parameters{<:Barotropic},
                                args...)
    return nothing
end

function initialize_orography!( orog::AbstractOrography,
                                ::Type{<:EarthOrography},
                                P::Parameters{M},
                                S::SpectralTransform,
                                G::Geometry) where {M<:Union{ShallowWater,PrimitiveEquation}}

    @unpack orography_path, orography_file, gravity = P
    @unpack Grid, lmax, mmax = S

    # LOAD NETCDF FILE (but not its data yet)
    if orography_path == ""
        path = joinpath(@__DIR__,"../input_data",orography_file)
    else
        path = joinpath(orography_path,orography_file)
    end
    ncfile = NetCDF.open(path)

    orography_highres = ncfile.vars["orog"][:,:]        # height [m]

    # Interpolate/coarsen to desired resolution
    #TODO also read lat,lon from file and flip array in case it's not as expected
    recompute_legendre = true
    orography_spec = spectral(orography_highres;Grid=FullGaussianGrid,recompute_legendre)
    
    @unpack orography, geopot_surf = orog
    copyto!(geopot_surf,orography_spec)     # truncates to the size of geopot_surf, no *gravity yet
    gridded!(orography,geopot_surf,S)       # to grid-point space
    geopot_surf .*= gravity                 # turn orography into surface geopotential
    spectral_truncation!(geopot_surf,lmax,mmax) # set the lmax+1 harmonics to zero

    # SCALE OROGRAPHY (or disable)
    orography .*= P.orography_scale
    geopot_surf .*= P.orography_scale
end