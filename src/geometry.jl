"""
    Geometry{NF<:AbstractFloat}

Geometry struct containing parameters and arrays describing an iso-latitude grid <:AbstractGrid
and the vertical levels. NF is the number format used for the precomputed constants.
"""
struct Geometry{NF<:AbstractFloat}      # NF: Number format

    # GRID TYPE AND RESOLUTION
    Grid::Type{<:AbstractGrid}
    nresolution::Int    # resolution parameter nlat_half or nside for HEALPix

    # GRID-POINT SPACE
    nlon::Int           # Maximum number of longitudes (at/around Equator)
    nlat::Int           # Number of latitudes
    nlev::Int           # Number of vertical levels
    nlat_half::Int      # Number of latitudes in one hemisphere (incl Equator)
    nside::Int          # HEALPix only, nside^2 are the # of grid points in each of the 12 base pixels
    npoints::Int        # total number of grid points
    radius_earth::Real  # Earth's radius [m]

    # LATITUDES (either Gaussian, equi-angle or HEALPix lats, depending on grid)
    lat::Vector{NF}         # array of latitudes (π/2...-π/2)
    latd::Vector{Float64}   # array of latitudes in degrees (90˚...-90˚)
    colat::Vector{NF}       # array of colatitudes (0...π)
    colatd::Vector{Float64} # array of colatitudes in degrees (0...180˚)

    # LONGITUDES
    lon::Vector{NF}         # full grids only: array of longitudes (0...2π)
    lond::Vector{Float64}   # array of longitudes in degrees (0...360˚)

    # COORDINATES
    lons::Vector{NF}        # longitude (0...2π) for each grid point in ring order
    lats::Vector{NF}        # latitude (π/2...-π/2) for each grid point in ring order

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    n_stratosphere_levels::Int      # number of upper levels for stratosphere
    σ_levels_half::Vector{NF}       # σ at half levels
    σ_levels_full::Vector{NF}       # σ at full levels
    σ_levels_thick::Vector{NF}      # σ level thicknesses
    σ_levels_half⁻¹_2::Vector{NF}   # 1/(2σ_levels_full)       #TODO rename?
    σ_f::Vector{NF}                 # akap/(2σ_levels_thick)   #TODO rename?

    # SINES AND COSINES OF LATITUDE
    sinlat::Vector{NF}              # sin of latitudes
    coslat::Vector{NF}              # cos of latitudes
    coslat⁻¹::Vector{NF}            # = 1/cos(lat)
    coslat²::Vector{NF}             # = cos²(lat)
    coslat⁻²::Vector{NF}            # = 1/cos²(lat)

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis::Vector{NF}          # = 2Ω*sin(lat)*radius_earth

    # GEOPOTENTIAL CALCULATION WORK ARRAYS
    xgeop1::Vector{NF}                  # ?
    xgeop2::Vector{NF}                  # ?
    lapserate_correction::Vector{NF}    # ?

    # TEMPORARY, development area. All these variables need to be checked for consistency and potentially defined somewhere else
    # tref ::Array{NF,1}   #temporarily defined here. Also defined in the Implict struct which is incomplete at the time of writing
    # rgas ::NF
    # fsgr ::Array{NF,1} 
    # tref3 ::Array{NF,1} 
end

"""
    G = Geometry(P::Parameters)

Generator function to create the Geometry struct from parameters in `P`.
"""
function Geometry(P::Parameters)

    @unpack Grid, trunc, nlev = P                   # grid type, spectral truncation, # of vertical levels
    @unpack radius_earth, rotation_earth, akap = P  # radius of earth, angular frequency, ratio of gas consts
    @unpack n_stratosphere_levels = P               # number of vertical levels used for stratosphere

    # RESOLUTION PARAMETERS
    nresolution = get_resolution(Grid,trunc)        # resolution parameter, nlat_half or nside for HEALPixGrid
    nlat_half = get_nlat_half(Grid,nresolution)     # contains equator for HEALPix
    nlat = 2nlat_half - nlat_odd(Grid)              # one less if grids have odd # of latitude rings
    nlon = get_nlon(Grid,nresolution)               # number of longitudes around the equator
    nside = Grid isa HEALPixGrid ? nresolution : 0  # nside is only defined for HEALPixGrid (npoints)
    npoints = get_npoints(Grid,nresolution)         # total number of grid points

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    colat = get_colat(Grid,nresolution)             # colatitude in radians
    lat = π/2 .- colat                              # latitude in radians
    colatd = colat*360/2π                           # and the same in degree 0...180˚
    latd = lat*360/2π                               # 90˚...-90˚

    # LONGITUDE VECTORS (empty for reduced grids)
    lon = get_lon(Grid,nresolution)                 # array of longitudes 0...2π
    lond = lon*360/2π                               # array of longitudes in degrees 0...360˚

    # COORDINATES for every grid point in ring order
    lats,lons = get_colatlons(Grid,nresolution)     # in radians

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
    coslat²  = coslat.^2    
    coslat⁻² = 1 ./ coslat²

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis = 2rotation_earth*sinlat*radius_earth

    # GEOPOTENTIAL coefficients to calculate geopotential (TODO reference)
    xgeop1 = zeros(nlev)
    xgeop2 = zeros(nlev)
    for k in 1:nlev
        xgeop1[k] = radius_earth*log(σ_levels_half[k+1]/σ_levels_half[k])
        if k != nlev
            xgeop2[k+1] = radius_earth*log(σ_levels_full[k+1]/σ_levels_half[k+1])
        end
    end

    if P.model == :primitive
        # LAPSE RATE correction (TODO reference)
        lapserate_correction = zeros(nlev-2)
        for k in 2:nlev-1
            lapserate_correction[k-1] = 0.5*xgeop1[k]*
                        log(σ_levels_half[k+1]/σ_levels_full[k]) / log(σ_levels_full[k+1]/σ_levels_full[k-1])
        end
    else
        lapserate_correction = zeros(1)
    end

    # Extra definitions. These will need to be defined consistently either here or somewhere else
    # Just defined here to proivide basic structure to allow for testing of other components of code
    # tref = 288.0max.(0.2, σ_levels_full) #more corrections needed here 
    # rgas = (2.0/7.0) / 1004.0
    # fsgr = (tref * 0.0) #arbitrary definition. Must be defined elsewhere 
    # tref3=fsgr.*tref #this actually is the correct definition. Needs better naming convention 

    # conversion to number format NF happens here
    Geometry{P.NF}( Grid,nresolution,
                    nlon,nlat,nlev,nlat_half,nside,npoints,radius_earth,
                    lat,latd,colat,colatd,lon,lond,lons,lats,
                    n_stratosphere_levels,
                    σ_levels_half,σ_levels_full,σ_levels_thick,σ_levels_half⁻¹_2,σ_f,
                    sinlat,coslat,coslat⁻¹,coslat²,coslat⁻²,
                    f_coriolis,
                    xgeop1,xgeop2,lapserate_correction)
                    # tref,rgas,fsgr,tref3)
end

"""
    σ_levels_half = vertical_coordinates(P::Parameters)

Vertical sigma coordinates defined by their nlev+1 half levels `σ_levels_half`. Sigma coordinates are
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