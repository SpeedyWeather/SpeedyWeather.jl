"""
    Geometry{NF<:AbstractFloat}

Geometry struct containing parameters and arrays describing an iso-latitude grid <:AbstractGrid
and the vertical levels. NF is the number format used for the precomputed constants.
"""
struct Geometry{NF<:AbstractFloat}      # NF: Number format

    # GRID TYPE AND RESOLUTION
    Grid::Type{<:AbstractGrid}
    nresolution::Int    # resolution parameter nlat_half

    # GRID-POINT SPACE
    nlon_max::Int       # Maximum number of longitudes (at/around Equator)
    nlon::Int           # Same (used for compatibility)
    nlat::Int           # Number of latitudes
    nlev::Int           # Number of vertical levels
    nlat_half::Int      # Number of latitudes in one hemisphere (incl Equator)
    npoints::Int        # total number of grid points
    radius_earth::Real  # Earth's radius [m]

    # LATITUDES (either Gaussian, equi-angle, HEALPix or HEALPix4 lats, depending on grid)
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

    # SINES AND COSINES OF LATITUDE
    sinlat::Vector{NF}              # sin of latitudes
    coslat::Vector{NF}              # cos of latitudes
    coslat⁻¹::Vector{NF}            # = 1/cos(lat)
    coslat²::Vector{NF}             # = cos²(lat)
    coslat⁻²::Vector{NF}            # = 1/cos²(lat)

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis::Vector{NF}          # = 2Ω*sin(lat)*radius_earth

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    n_stratosphere_levels::Int      # number of upper levels for stratosphere
    σ_levels_half::Vector{NF}       # σ at half levels
    σ_levels_full::Vector{NF}       # σ at full levels
    σ_levels_thick::Vector{NF}      # σ level thicknesses
    σ_levels_thick⁻¹_half::Vector{NF}   # = 1/(2σ_levels_thick)
    σ_f::Vector{NF}                 # akap/(2σ_levels_thick)   #TODO rename?
    σ_lnp_A::Vector{NF}
    σ_lnp_B::Vector{NF}

    # VERTICAL REFERENCE TEMPERATURE PROFILE
    temp_ref_profile::Vector{NF}
    # temp_ref_profile_σ::Vector{NF}

    # GEOPOTENTIAL INTEGRATION (on half/full levels)
    Δp_geopot_half::Vector{NF}      # = R*(ln(p_k+1) - ln(p_k+1/2)), for half level geopotential
    Δp_geopot_full::Vector{NF}      # = R*(ln(p_k+1/2) - ln(p_k)), for full level geopotential
    lapserate_corr::Vector{NF}      # ?

    # PARAMETERIZATIONS
    entrainment_profile::Vector{NF}
end

"""
    G = Geometry(P::Parameters)

Generator function to create the Geometry struct from parameters in `P`.
"""
function Geometry(P::Parameters,Grid::Type{<:AbstractGrid})

    @unpack trunc, nlev = P                         # grid type, spectral truncation, # of vertical levels
    @unpack radius_earth, rotation_earth = P        # radius of earth, angular frequency
    @unpack R_dry, cₚ = P                           # gas constant for dry air, heat capacity
    @unpack n_stratosphere_levels = P               # number of vertical levels used for stratosphere
    @unpack temp_ref, temp_top, lapse_rate, gravity = P # for reference atmosphere

    # RESOLUTION PARAMETERS
    nresolution = get_resolution(Grid,trunc)        # resolution parameter nlat_half
    nlat_half = nresolution                         # number of latitude rings on one hemisphere (Equator incl)
    nlat = get_nlat(Grid,nlat_half)                 # 2nlat_half but one less if grids have odd # of lat rings
    nlon_max = get_nlon_max(Grid,nlat_half)         # number of longitudes around the equator
    nlon = nlon_max                                 # same (used for compatibility)
    npoints = get_npoints(Grid,nlat_half)           # total number of grid points

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    colat = get_colat(Grid,nlat_half)               # colatitude in radians
    lat = π/2 .- colat                              # latitude in radians
    colatd = colat*360/2π                           # and the same in degree 0...180˚
    latd = lat*360/2π                               # 90˚...-90˚

    # LONGITUDE VECTORS (empty for reduced grids)
    lon = get_lon(Grid,nlat_half)                   # array of longitudes 0...2π (full grids only)
    lond = lon*360/2π                               # array of longitudes in degrees 0...360˚

    # COORDINATES for every grid point in ring order
    lats,lons = get_colatlons(Grid,nlat_half)       # in radians

    # SINES AND COSINES OF LATITUDE
    sinlat = sin.(lat)
    coslat = cos.(lat)
    coslat⁻¹ = 1 ./ coslat
    coslat²  = coslat.^2
    coslat⁻² = 1 ./ coslat²

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis = 2rotation_earth*sinlat*radius_earth

    # VERTICAL SIGMA COORDINATE
    # σ = p/p0 (fraction of surface pressure)
    # sorted such that σ_levels_half[end] is at the planetary boundary
    κ = R_dry/cₚ
    σ_levels_half = vertical_coordinates(P)
    σ_levels_full = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])
    σ_levels_thick = σ_levels_half[2:end] - σ_levels_half[1:end-1]
    σ_levels_thick⁻¹_half = 1 ./ (2σ_levels_thick)
    σ_f = κ ./ (2σ_levels_full)

    σ_lnp_A = [log(σ_levels_full[k]/σ_levels_half[k])/σ_levels_thick[k] for k in 1:nlev]
    σ_lnp_B = [log(σ_levels_half[k+1]/σ_levels_full[k])/σ_levels_thick[k] for k in 1:nlev]

    # VERTICAL REFERENCE TEMPERATURE PROFILE
    RΓ = R_dry*lapse_rate/(1000*gravity)                    # origin unclear but profile okay
    temp_ref_profile = [max(temp_top,temp_ref*σ^RΓ) for σ in σ_levels_full]
    # temp_ref_profile_σ = temp_ref_profile .* σ_f

    # GEOPOTENTIAL, coefficients to calculate geopotential
    Δp_geopot_half, Δp_geopot_full = initialise_geopotential(σ_levels_full,σ_levels_half,R_dry)

    # LAPSE RATE correction
    lapserate_corr = lapserate_correction(σ_levels_full,σ_levels_half,Δp_geopot_full)

    # Compute the entrainment coefficients for the convection parameterization.
    @unpack max_entrainment = P
    entrainment_profile = zeros(nlev)
    for k = 2:nlev-1
        entrainment_profile[k] = max(0, (σ_levels_full[k] - 0.5)^2)
    end

    # profile as fraction of cloud-base mass flux
    entrainment_profile /= sum(entrainment_profile)  # Normalise
    entrainment_profile *= max_entrainment           # fraction of max entrainment

    # conversion to number format NF happens here
    Geometry{P.NF}( Grid,nresolution,
                    nlon_max,nlon,nlat,nlev,nlat_half,npoints,radius_earth,
                    lat,latd,colat,colatd,lon,lond,lons,lats,
                    sinlat,coslat,coslat⁻¹,coslat²,coslat⁻²,f_coriolis,
                    n_stratosphere_levels,
                    σ_levels_half,σ_levels_full,σ_levels_thick,σ_levels_thick⁻¹_half,σ_f,
                    σ_lnp_A,σ_lnp_B,
                    temp_ref_profile,
                    Δp_geopot_half,Δp_geopot_full,lapserate_corr,entrainment_profile)
                    # tref,rgas,fsgr,tref3)
end

# use Grid in Parameters if not provided
Geometry(P::Parameters) = Geometry(P,P.Grid)

"""
    σ_levels_half = vertical_coordinates(P::Parameters)

Vertical sigma coordinates defined by their nlev+1 half levels `σ_levels_half`. Sigma coordinates are
fraction of surface pressure (p/p0) and are sorted from top (stratosphere) to bottom (surface).
The first half level is at 0 the last at 1. Evaluate a generalised logistic function with
coefficients in `P` for the distribution of values in between. Default coefficients follow
the L31 configuration historically used at ECMWF."""
function vertical_coordinates(P::Parameters)
    @unpack nlev,GLcoefs,σ_levels_half = P

    if length(σ_levels_half) == 0       # choose σ levels automatically
        halflevels_normalised = range(0,1,nlev+1)   # normalised = level/nlev
        σ_levels_half = generalised_logistic(halflevels_normalised,GLcoefs)
        σ_levels_half[1] = 0            # topmost half-level is at 0 pressure
        σ_levels_half[end] = 1          # lowermost half-level is at p=p_surface
    else                                # choose σ levels manually
        @assert σ_levels_half[1] == 0 "First manually specified σ_levels_half has to be zero 0"
        @assert σ_levels_half[end] == 1 "Last manually specified σ_levels_half has to be 1."
        @assert nlev == (length(σ_levels_half) - 1) "nlev has to be length of σ_levels_half - 1"
    end

    @assert isincreasing(σ_levels_half) "Vertical coordinates are not increasing."
    return σ_levels_half
end

"""Generalised logistic function based on the coefficients in `coefs`."""
function generalised_logistic(x,coefs::GenLogisticCoefs)
    @unpack A,K,C,Q,B,M,ν = coefs
    return @. A + (K-A)/(C+Q*exp(-B*(x-M)))^inv(ν)
end
