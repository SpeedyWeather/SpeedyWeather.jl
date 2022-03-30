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
    lon::Array{NF,1}    # array of longitudes
    lat::Array{NF,1}    # array of latitudes

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
    sinlat_NH::Array{NF,1}      # only northern hemisphere
    coslat_NH::Array{NF,1}      # only northern hemisphere
    radang::Array{NF,1}         # radians of latitudes TODO rename to radlat?
    cosg::Array{NF,1}           # this should be sinlat?!
    cosg⁻¹::Array{NF,1}         # rename to sinlat⁻¹?
    cosg⁻²::Array{NF,1}         # rename to sinlat⁻²?

    # CORIOLIS FREQUENCY
    f_coriolis::Array{NF,1}              # = 2Ω*sin(lat)

    # GEOPOTENTIAL CALCULATION WORK ARRAYS
    xgeop1::Array{NF,1}                  # ?
    xgeop2::Array{NF,1}                  # ?
    lapserate_correction::Array{NF,1}    # ?


    #TEMPORARY, development area. All these variables need to be checked for consistency and potentially defined somewhere else
    tref ::Array{NF,1}   #temporarily defined here. Also defined in the Implict struct which is incomplete at the time of writing
    rgas ::NF
    fsgr ::Array{NF,1} 
    tref3 ::Array{NF,1} 



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

    # GRID SPACE ARRAYS GAUSSIAN GRID - lon is equi-spaced, lat is not!
    dlon = 360 / nlon                       # grid spacing in longitude
    dlat = 180 / nlat                       # average grid spacing in latitude
    lon  = Array(0:dlon:360-dlon)           # array of longitudes
    # array of latitudes (North to South) corresponding to the zeros
    # of the (unassociated) legendre polynomial order nlat
    lat  = reverse(asind.(FastGaussQuadrature.gausslegendre(nlat)[1]))

    # VERTICAL SIGMA COORDINATE 
    # σ = p/p0 (fraction of surface pressure)
    # sorted such that σ_levels_half[end] is at the planetary boundary
    σ_levels_half = vertical_coordinates(P)
    σ_levels_full = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])
    σ_levels_thick = σ_levels_half[2:end] - σ_levels_half[1:end-1]
    σ_levels_half⁻¹_2 = 1 ./ (2σ_levels_thick)
    σ_f = akap ./ (2σ_levels_full)


    # SINES AND COSINES OF LATITUDE
    sinlat = sind.(lat)
    coslat = cosd.(lat)
    sinlat_NH = sinlat[1:nlat_half] # sinlat only for northern hemisphere = NH
    coslat_NH = coslat[1:nlat_half]
    radang = asin.(sinlat)
    cosg   = sinlat                 # inconsistent here due to the sin/cos swap
    cosg⁻¹ = 1 ./ cosg
    cosg⁻² = 1 ./ cosg.^2

    # CORIOLIS FREQUENCY
    f_coriolis = 2Ω*sinlat

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


    #Extra definitions. These will need to be defined consistently either here or somewhere else
    #Just defined here to proivide basic structure to allow for testing of other components of code


    tref = 288.0max.(0.2, σ_levels_full) #more corrections needed here 
    rgas = (2.0/7.0) / 1004.0
    fsgr = (tref * 0.0) #arbitrary definition. Must be defined elsewhere 
    tref3=fsgr.*tref #this actually is the correct definition. Needs better naming convention 


    # conversion to number format NF happens here
    Geometry{P.NF}( nlon,nlat,nlev,nlat_half,nlon_half,
                    dlon,dlat,lon,lat,
                    n_stratosphere_levels,
                    σ_levels_half,σ_levels_full,σ_levels_thick,σ_levels_half⁻¹_2,σ_f,
                    sinlat,coslat,sinlat_NH,coslat_NH,radang,
                    cosg,cosg⁻¹,cosg⁻²,f_coriolis,xgeop1,xgeop2,lapserate_correction,
                    tref,rgas,fsgr,tref3)
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