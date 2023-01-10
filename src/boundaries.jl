# gridded field (orography) and spectral field (surface geopotential) for Earth's orography
struct EarthOrography{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF}
    orography::Grid                                 # height [m]
    geopot_surf::LowerTriangularMatrix{Complex{NF}} # surface geopotential, height*gravity [m²/s²]
end

abstract type ZonalRidge{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF} end
ZonalRidge(orography,geopot_surf) = EarthOrography(orography,geopot_surf)

abstract type NoOrography{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractOrography{NF} end
NoOrography(orography,geopot_surf) = EarthOrography(orography,geopot_surf)

function Base.zeros(O::Type{<:AbstractOrography},S::SpectralTransform{NF}) where NF
    @unpack Grid, nlat_half, lmax, mmax = S
    orography   = zeros(Grid{NF},nlat_half)
    geopot_surf = zeros(LowerTriangularMatrix{Complex{NF}},lmax+2,mmax+1)
    return O(orography,geopot_surf)
end

struct Boundaries{NF<:AbstractFloat,Grid<:AbstractGrid{NF}} <: AbstractBoundaries{NF}
    orography::AbstractOrography{NF}
    # landsea_mask::AbstractLandSeaMask{NF}
    # albedo::AbstractAlbedo{NF}
end

function Boundaries(P::Parameters,
                    S::SpectralTransform{NF},
                    G::Geometry{NF}) where NF

    orography = zeros(P.orography,S)
    initialize_orography!(orography,P.orography,P,S,G)
    return Boundaries{NF,S.Grid{NF}}(orography)
end

# recalculate spectral transform or geometry if not provided
Boundaries(P::Parameters,S::SpectralTransform) = Boundaries(P,S,Geometry(P))
Boundaries(P::Parameters) = Boundaries(P,SpectralTransform(P))

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
    @unpack lmax, mmax = S

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
    recompute_legendre = true   # don't allocate large arrays as spectral transform is not reused
    Grid = FullGaussianGrid     # grid the orography file comes with
    orography_spec = spectral(orography_highres;Grid,recompute_legendre)
    
    @unpack orography, geopot_surf = orog
    copyto!(geopot_surf,orography_spec)     # truncates to the size of geopot_surf, no *gravity yet
    gridded!(orography,geopot_surf,S)       # to grid-point space
    geopot_surf .*= gravity                 # turn orography into surface geopotential
    spectral_truncation!(geopot_surf,lmax,mmax) # set the lmax+1 harmonics to zero

    # SCALE OROGRAPHY (or disable)
    orography .*= P.orography_scale
    geopot_surf .*= P.orography_scale
end