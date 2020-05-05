struct Geometry{T<:AbstractFloat}

    # GRID-POINT SPACE
    nlon::Int           # Number of longitudes
    nlat::Int           # Number of latitudes
    nlev::Int           # Number of vertical levels
    nlat_half::Int      # Number of latitudes in one hemisphere

    dlon::T             # grid spacing in longitude
    dlat::T             # average grid spacing in latitude
    lon::Array{T,1}     # array of longitudes
    lat::Array{T,1}     # array of latitudes

    # VERTICAL SIGMA COORDINATE σ = p/p₀ (fraction of surface pressure)
    σ_half::Array{T,1}      # σ at half levels
    σ_full::Array{T,1}      # σ at full levels
    σ_thick::Array{T,1}     # σ thicknesses
    σ_half⁻¹_2::Array{T,1}  # 1/(2σ_full)       #TODO rename?
    σ_f::Array{T,1}         # akap/(2σ_thick)   #TODO rename?

    # SINES AND COSINES OF LATITUDE
    sinlat::Array{T,1}          # sin of latitudes
    coslat::Array{T,1}          # cos of latitudes
    sinlat_NH::Array{T,1}       # only northern hemisphere
    coslat_NH::Array{T,1}       # only northern hemisphere
    radang::Array{T,1}          # radians of latitudes TODO rename to radlat?
    cosg::Array{T,1}            # this should be sinlat?!
    cosg⁻¹::Array{T,1}          # rename to sinlat⁻¹?
    cosg⁻²::Array{T,1}          # rename to sinlat⁻²?

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
function Geometry{T}(C::Constants,P::Params) where T

    @unpack R,Ω,akap = C
    @unpack nlon,nlat,nlev,trunc = P

    nlat_half = nlat ÷ 2

    # GRID SPACE ARRAYS GAUSSIAN GRID - lon is equi-spaced, lat is not!
    dlon = 360 / nlon                       # grid spacing in longitude
    dlat = 180 / nlat                       # average grid spacing in latitude
    lon  = Array(0:dlon:360-dlon)           # array of longitudes
    lat  = asind.(gausslegendre(nlat)[1])   # array of latitudes
                                            # which correspond to the zeros
                                            # of the legendre polynomial order nlat

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    # sorted such that σ_half[end] is at the planetary boundary
    #TODO make nlev-dependent
    σ_half = [0.0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0]
    σ_full = 0.5*(σ_half[2:end] + σ_half[1:end-1])
    σ_thick = σ_half[2:end] - σ_half[1:end-1]
    σ_half⁻¹_2 = 1 ./ (2σ_thick)
    σ_f = akap ./ (2σ_full)

    # SINES AND COSINES OF LATITUDE
    sinlat = sind.(lat)
    coslat = cosd.(lat)
    sinlat_NH = sinlat[nlat_half+1:end]     # sinlat only for northern hemisphere
    coslat_NH = coslat[nlat_half+1:end]
    radang = asin.(sinlat)
    cosg   = sinlat                         # inconsistent here due to the sin/cos swap
    cosg⁻¹ = T(1.0)./cosg
    cosg⁻² = T(1.0)./cosg.^T(2.0)

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

    # LAPSE RATE
    lapserate_correction = zeros(nlev-2)
    for k in 2:nlev-1
        lapserate_correction[k-1] = 0.5*xgeop1[k]*
                    log(σ_half[k+1]/σ_full[k]) / log(σ_full[k+1]/σ_full[k-1])
    end

    # conversion to T happens here
    Geometry{T}(nlon,nlat,nlev,nlat_half,dlon,dlat,lon,lat,
                σ_half,σ_full,σ_thick,σ_half⁻¹_2,σ_f,
                sinlat,coslat,sinlat_NH,coslat_NH,radang,
                cosg,cosg⁻¹,cosg⁻²,f,xgeop1,xgeop2,lapserate_correction)
end
