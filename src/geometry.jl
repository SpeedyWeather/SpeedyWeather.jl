"""
Geometry struct containing parameters and arrays describing the Gaussian grid
and the vertical levels. NF is the number format used for the precomputed constants.
"""
struct Geometry{NF<:AbstractFloat}      # NF: Number format

    # GRID-POINT SPACE
    nlon::Int           # Number of longitudes
    nlat::Int           # Number of latitudes
    nlev::Int           # Number of vertical levels
    nlat_half::Int      # Number of latitudes in one hemisphere
    nlon_half::Int      # Half the number of longitudes

    dlon::NF            # grid spacing in longitude
    dlat::NF            # average grid spacing in latitude
    lon::Array{NF,1}    # array of longitudes (0...2π)
    lond::Array{NF,1}   # array of longitudes in degrees (0...360˚)
    lat::Array{NF,1}    # array of latitudes (π/2...-π/2)
    latd::Array{NF,1}   # array of latitudes in degrees (90˚...-90˚)
    colat::Array{NF,1}  # array of colatitudes (0...π)
    colatd::Array{NF,1} # array of colatitudes in degrees (0...180˚)

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    n_stratosphere_levels::Int      # number of upper levels for stratosphere
    σ_levels_half::Array{NF,1}      # σ at half levels
    σ_levels_full::Array{NF,1}      # σ at full levels
    σ_levels_thick::Array{NF,1}     # σ level thicknesses
    σ_levels_half⁻¹_2::Array{NF,1}  # 1/(2σ_levels_full)       #TODO rename?
    σ_f::Array{NF,1}                # akap/(2σ_levels_thick)   #TODO rename?

    # SINES AND COSINES OF LATITUDE
    sinlat::Array{NF,1}         # sin of latitudes
    coslat::Array{NF,1}         # cos of latitudes
    coslat⁻¹::Array{NF,1}       # =1/cos(lat)

    # CORIOLIS FREQUENCY
    f::Array{NF,1}              # = 2Ω*sin(lat)

    # GEOPOTENTIAL CALCULATION WORK ARRAYS
    xgeop1::Array{NF,1}                  # ?
    xgeop2::Array{NF,1}                  # ?
    lapserate_correction::Array{NF,1}    # ?
end

"""
Defines the geometry.
"""
function Geometry(P::Parameters)

    @unpack nlon,nlat,nlev,trunc = P
    @unpack R,Ω,akap = P
    @unpack n_stratosphere_levels = P

    nlat_half = nlat ÷ 2
    nlon_half = nlon ÷ 2

    # GRID SPACE ARRAYS lon is equi-spaced
    dlon = 360 / nlon                       # grid spacing in longitude
    dlat = 180 / nlat                       # average grid spacing in latitude
    lond = Array(0:dlon:360-dlon)           # array of longitudes in degrees 0...360˚
    lon = lond/360*2π                       # array of longitudes 0...2π
    
    # REGULAR GAUSSIAN GRID, latitudes are not equi-distant
    # Zero nodes of the (unassociated) legendre polynomial order nlat
    nodes = FastGaussQuadrature.gausslegendre(nlat)[1]
    colat = π .- acos.(nodes)               # colatitudes 0...π
    colatd = 180 .- cosd.(nodes)            # colatitudes in degrees 0...180˚
    lat = -sin.(nodes)                      # latitudes π/2...-π/2
    latd = -sind.(nodes)                    # latitudes in degrees 90˚...-90˚

    # VERTICAL SIGMA COORDINATE 
    # σ = p/p0 (fraction of surface pressure)
    # sorted such that σ_levels_half[end] is at the planetary boundary
    σ_levels_half = vertical_coordinates(P)
    σ_levels_full = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])
    σ_levels_thick = σ_levels_half[2:end] - σ_levels_half[1:end-1]
    σ_levels_half⁻¹_2 = 1 ./ (2σ_levels_thick)
    σ_f = akap ./ (2σ_levels_full)

    # SINES AND COSINES OF LATITUDE
    sinlat = sin.(lat)
    coslat = cos.(lat)
    coslat⁻¹ = 1 ./ coslat

    # CORIOLIS FREQUENCY
    f = 2Ω*sinlat

    # GEOPOTENTIAL coefficients to calculate geopotential (TODO reference)
    xgeop1 = zeros(nlev)
    xgeop2 = zeros(nlev)
    for k in 1:nlev
        xgeop1[k] = R*log(σ_levels_half[k+1]/σ_levels_half[k])
        if k != nlev
            xgeop2[k+1] = R*log(σ_levels_full[k+1]/σ_levels_half[k+1])
        end
    end

    # LAPSE RATE correction (TODO reference)
    lapserate_correction = zeros(nlev-2)
    for k in 2:nlev-1
        lapserate_correction[k-1] = 0.5*xgeop1[k]*
                    log(σ_levels_half[k+1]/σ_levels_full[k]) / log(σ_levels_full[k+1]/σ_levels_full[k-1])
    end

    # conversion to number format NF happens here
    Geometry{P.NF}( nlon,nlat,nlev,nlat_half,nlon_half,
                    dlon,dlat,lon,lond,lat,latd,colat,colatd,
                    n_stratosphere_levels,
                    σ_levels_half,σ_levels_full,σ_levels_thick,σ_levels_half⁻¹_2,σ_f,
                    sinlat,coslat,coslat⁻¹,
                    f,
                    xgeop1,xgeop2,lapserate_correction)
end

"""Vertical sigma coordinates defined by their nlev+1 half levels `σ_levels_half`. Sigma coordinates are
fraction of surface pressure (p/p0) and are sorted from top (stratosphere) to bottom (surface).
The first half level is at 0 the last at 1. Evaluate a generalised logistic function with
coefficients in `P` for the distribution of values in between. Default coefficients follow
the L31 configuration historically used at ECMWF."""
function vertical_coordinates(P::Parameters)
    @unpack nlev,GLcoefs = P

    halflevels_normalised = range(0,1,nlev+1)   # normalised = level/nlev 
    σ_levels_half = generalised_logistic(halflevels_normalised,GLcoefs)
    σ_levels_half[1] = 0           # topmost half-level is at 0 pressure
    σ_levels_half[end] = 1         # lowermost half-level is at p=p_surface

    @assert isincreasing(σ_levels_half) "Vertical coordinates are not increasing."
    return σ_levels_half
end

"""Generalised logistic function based on the coefficients in `coefs`."""
function generalised_logistic(x,coefs::GenLogisticCoefs)
    @unpack A,K,C,Q,B,M,ν = coefs
    return @. A + (K-A)/(C+Q*exp(-B*(x-M)))^(1/ν)
end