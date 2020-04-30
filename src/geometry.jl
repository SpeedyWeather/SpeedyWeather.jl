struct Geometry{T<:AbstractFloat}

    # GRID-POINT SPACE
    nlon::Int           # Number of longitudes
    nlat::Int           # Number of latitudes
    nlev::Int           # Number of vertical levels

    dlat::T             # grid spacing in latitude
    dlon::T             # grid spacing in longitude
    lat::Array{T,1}     # array of latitudes
    lon::Array{T,1}     # array of longitudes

    # SPECTRAL SPACE
    trunc::Int      # Spectral truncation
    nx::Int         # Number of total wavenumbers
    mx::Int         # Number of zonal wavenumbers

    # VERTICAL SIGMA COORDINATE σ = p/p₀ (fraction of surface pressure)
    σ_half::Array{T,1}      # σ at half levels
    σ_full::Array{T,1}      # σ at full levels
    σ_thick::Array{T,1}     # σ thicknesses
    σ_half⁻¹_2::Array{T,1}  # 1/(2σ_full)       #TODO rename?
    σ_f::Array{T,1}         # akap/(2σ_thick)   #TODO rename?

    # SINES AND COSINES OF LATITUDE
    sinlat::Array{T,1}
    coslat::Array{T,1}
    sinlat_half::Array{T,1}
    coslat_half::Array{T,1}
    radang::Array{T,1}
    cosg::Array{T,1}            # ?
    cosg⁻¹::Array{T,1}          # ?
    cosg⁻²::Array{T,1}          # ?

    # CORIOLIS FREQUENCY
    f::Array{T,1}               # = 2Ω*sin(lat)

    # GEOPOTENTIAL CALCULATION WORK ARRAYS
    xgeop1::Array{T,1}                  # ?
    xgeop2::Array{T,1}                  # ?
    lapserate_correction::Array{T,1}    # ?
end

"""
Defines the geometry.
"""
function Geometry{T}(   nlon::Int,
                        nlat::Int,
                        nlev::Int,
                        trunc::Int,
                        constants::Constants) where T

    @unpack R,Ω,akap = constants

    @assert nlev == 8   "Only 8 model levels are currently supported"
    @assert nlat == 48  "Only nlat=48 currently supported"
    @assert nlon == 96  "Only nlon=96 currently supported"

    # GRID SPACE ARRAYS
    dlat = 180 / nlat                   # grid spacing in latitude
    dlon = 360 / nlon                   # grid spacing in longitude
    lat = Array(-90+dlat/2:dlat:90-dlat/2)    # array of latitudes
    lon = Array(0:dlon:360-dlon)              # array of longitudes

    # VERTICAL SIGMA COORDINATE σ = p/p₀ (fraction of surface pressure)
    # sorted such that σ_half[end] is at the planetary boundary
    #TODO make nlev-dependent
    σ_half = [0.0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0]
    σ_full = 0.5*(σ_half[2:end] + σ_half[1:end-1])
    σ_thick = σ_half[2:end] - σ_half[1:end-1]
    σ_half⁻¹_2 = 1 ./ (2σ_thick)
    σ_f = akap ./ (2.0σ_full)

    # SPECTRAL
    # TODO make trunc a function of nlat,nlon?
    nx = trunc+2
    mx = trunc+1

    # SINES AND COSINES OF LATITUDE
    sinlat = sind.(lat)
    coslat = cosd.(lat)
    sinlat_half = sinlat[nlat÷2+1:end]      # only Northern Hemisphere #TODO rename to sinlat_NH?
    coslat_half = coslat[nlat÷2+1:end]      # only Northern Hemisphere
    radang = 2π*lat/360                     # latitude in radians
    cosg   = coslat                         # ?
    cosg⁻¹ = 1 ./ cosg                      # ?
    cosg⁻² = cosg⁻¹.^2                      # ?

    # CORIOLIS FREQUENCY
    f = 2Ω*sinlat

    # GEOPOTENTIAL
    xgeop1 = zeros(nlev)                    # coefficients to calculate geopotential
    xgeop2 = zeros(nlev)                    # coefficients to calculate geopotential
    for k in 1:nlev
        xgeop1[k] = R*log(σ_half[k+1]/σ_half[k])
        if k != nlev
            xgeop2[k+1] = R*log(σ_full[k+1]/σ_half[k+1])
        end
    end

    lapserate_correction = zeros(nlev-2)
    for k in 2:nlev-1
        lapserate_correction[k-1] = 0.5*xgeop1[k]*
                    log(σ_half[k+1]/σ_full[k]) / log(σ_full[k+1]/σ_full[k-1])
    end

    # conversion to T happens here
    Geometry{T}(nlon,nlat,nlev,dlat,dlon,lat,lon,
                trunc,nx,mx,
                σ_half,σ_full,σ_thick,σ_half⁻¹_2,σ_f,
                sinlat,coslat,sinlat_half,coslat_half,radang,
                cosg,cosg⁻¹,cosg⁻²,f,xgeop1,xgeop2,lapserate_correction)
end
