"""
    Geometry{NF<:AbstractFloat} <: AbstractGeometry

Geometry struct containing parameters and arrays describing an iso-latitude grid <:AbstractGrid
and the vertical levels. NF is the number format used for the precomputed constants.
"""
struct Geometry{NF<:AbstractFloat} <: AbstractGeometry{NF}     # NF: Number format

    # GRID TYPE AND RESOLUTION
    Grid::Type{<:AbstractGrid}
    nlat_half::Int      # resolution parameter nlat_half, # of latitudes on one hemisphere (incl Equator)

    # GRID-POINT SPACE
    nlon_max::Int       # Maximum number of longitudes (at/around Equator)
    nlon::Int           # Same (used for compatibility)
    nlat::Int           # Number of latitudes
    nlev::Int           # Number of vertical levels
    npoints::Int        # total number of grid points
    radius::Float64     # Planet's radius [m]

    # LATITUDES (either Gaussian, equi-angle, HEALPix or OctaHEALPix lats, depending on grid)
    latd::Vector{Float64}   # array of latitudes in degrees (90˚...-90˚)

    # LONGITUDES
    lond::Vector{Float64}   # array of longitudes in degrees (0...360˚), empty for non-full grids

    # COORDINATES
    londs::Vector{NF}       # longitude (-180˚...180˚) for each grid point in ring order
    latds::Vector{NF}       # latitude (-90˚...˚90) for each grid point in ring order

    # SINES AND COSINES OF LATITUDE
    sinlat::Vector{NF}              # sin of latitudes
    coslat::Vector{NF}              # cos of latitudes
    coslat⁻¹::Vector{NF}            # = 1/cos(lat)
    coslat²::Vector{NF}             # = cos²(lat)
    coslat⁻²::Vector{NF}            # = 1/cos²(lat)

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis::Vector{NF}          # = 2Ω*sin(lat)*radius

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    n_stratosphere_levels::Int      # number of upper levels for stratosphere
    σ_levels_half::Vector{NF}       # σ at half levels
    σ_levels_full::Vector{NF}       # σ at full levels
    σ_levels_thick::Vector{NF}      # σ level thicknesses
    ln_σ_levels_full::Vector{NF}    # log of σ at full levels

    σ_lnp_A::Vector{NF}             # σ-related factors needed for adiabatic terms 1st
    σ_lnp_B::Vector{NF}             # and 2nd

    # VERTICAL REFERENCE TEMPERATURE PROFILE
    temp_ref_profile::Vector{NF}

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

    (;trunc, dealiasing, nlev) = P             # grid type, spectral truncation, # of vertical levels
    (;radius, rotation, gravity) = P.planet    # radius of planet, angular frequency, gravity
    (;R_dry, cₚ) = P                           # gas constant for dry air, heat capacity
    (;σ_tropopause) = P                        # number of vertical levels used for stratosphere
    (;temp_ref, temp_top, lapse_rate) = P      # for reference atmosphere
    (;ΔT_stratosphere) = P                     # used for stratospheric temperature increase

    # RESOLUTION PARAMETERS
    # resolution parameter nlat_half (= number of lat rings on one hemisphere (Equator incl) 
    # from spectral resolution and dealiasing parameter (e.g. quadratic grid for T31)
    nlat_half = SpeedyTransforms.get_nlat_half(trunc,dealiasing)
    nlat = get_nlat(Grid,nlat_half)                 # 2nlat_half but -1 if grids have odd # of lat rings
    nlon_max = get_nlon_max(Grid,nlat_half)         # number of longitudes around the equator
    nlon = nlon_max                                 # same (used for compatibility)
    npoints = get_npoints(Grid,nlat_half)           # total number of grid points

    # LATITUDE VECTORS (based on Gaussian, equi-angle or HEALPix latitudes)
    latd = get_latd(Grid,nlat_half)                 # latitude in 90˚...-90˚

    # LONGITUDE VECTORS (empty for reduced grids)
    lond = get_lond(Grid,nlat_half)                 # array of longitudes 0...360˚ (full grids only)

    # COORDINATES for every grid point in ring order
    latds,londs = get_latdlonds(Grid,nlat_half)     # in -90˚...90˚N, -180˚...180˚

    # SINES AND COSINES OF LATITUDE
    sinlat = sind.(latd)
    coslat = cosd.(latd)
    coslat⁻¹ = 1 ./ coslat
    coslat²  = coslat.^2
    coslat⁻² = 1 ./ coslat²

    # CORIOLIS FREQUENCY (scaled by radius as is vorticity)
    f_coriolis = 2rotation*sinlat*radius

    # VERTICAL SIGMA COORDINATE
    # σ = p/p0 (fraction of surface pressure)
    # sorted such that σ_levels_half[end] is at the planetary boundary
    σ_levels_half = vertical_coordinates(P)
    σ_levels_full = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])
    σ_levels_thick = σ_levels_half[2:end] - σ_levels_half[1:end-1]
    ln_σ_levels_full = log.(vcat(σ_levels_full,1))              # include surface (σ=1) as last element

    # ADIABATIC TERM, Simmons and Burridge, 1981, eq. 3.12
    # precompute ln(σ_k+1/2) - ln(σ_k-1/2) but swap sign, include 1/Δσₖ
    σ_lnp_A = log.(σ_levels_half[1:end-1]./σ_levels_half[2:end]) ./ σ_levels_thick
    σ_lnp_A[1] = 0  # the corresponding sum is 1:k-1 so 0 to replace log(0) from above
    
    # precompute the αₖ = 1 - p_k-1/2/Δpₖ*log(p_k+1/2/p_k-1/2) term in σ coordinates
    σ_lnp_B = 1 .- σ_levels_half[1:end-1]./σ_levels_thick .*
                    log.(σ_levels_half[2:end]./σ_levels_half[1:end-1])
    σ_lnp_B[1] = σ_levels_half[1] <= 0 ? log(2) : σ_lnp_B[1]    # set α₁ = log(2), eq. 3.19
    σ_lnp_B .*= -1  # absorb sign from -1/Δσₖ only, eq. 3.12

    # TROPOPAUSE/STRATOSPHERIC LEVELS
    n_stratosphere_levels = sum(σ_levels_full .< σ_tropopause)  # of levels above σ_tropopause

    # VERTICAL REFERENCE TEMPERATURE PROFILE
    # integrate hydrostatic equation from pₛ to p, use ideal gas law p = ρRT and linear
    # temperature decrease with height: T = Tₛ - ΔzΓ with lapse rate Γ
    # for stratosphere (σ < σ_tropopause) increase temperature (Jablonowski & Williamson. 2006, eq. 5)
    RΓg⁻¹ = R_dry*lapse_rate/(1000*gravity)     # convert lapse rate from [K/km] to [K/m]
    temp_ref_profile = [temp_ref*σ^RΓg⁻¹ for σ in σ_levels_full]
    temp_ref_profile .+= [σ < σ_tropopause ? ΔT_stratosphere*(σ_tropopause-σ)^5 : 0 for σ in σ_levels_full]

    # GEOPOTENTIAL, coefficients to calculate geopotential
    Δp_geopot_half, Δp_geopot_full = initialize_geopotential(σ_levels_full,σ_levels_half,R_dry)

    # LAPSE RATE correction
    lapserate_corr = lapserate_correction(σ_levels_full,σ_levels_half,Δp_geopot_full)

    # Compute the entrainment coefficients for the convection parameterization.
    (;max_entrainment) = P
    entrainment_profile = zeros(nlev)
    for k = 2:nlev-1
        entrainment_profile[k] = max(0, (σ_levels_full[k] - 0.5)^2)
    end

    # profile as fraction of cloud-base mass flux
    entrainment_profile /= sum(entrainment_profile)  # Normalise
    entrainment_profile *= max_entrainment           # fraction of max entrainment

    # conversion to number format NF happens here
    Geometry{P.NF}( Grid,nlat_half,
                    nlon_max,nlon,nlat,nlev,npoints,radius,
                    latd,lond,londs,latds,
                    sinlat,coslat,coslat⁻¹,coslat²,coslat⁻²,f_coriolis,
                    n_stratosphere_levels,
                    σ_levels_half,σ_levels_full,σ_levels_thick,ln_σ_levels_full,
                    σ_lnp_A,σ_lnp_B,
                    temp_ref_profile,
                    Δp_geopot_half,Δp_geopot_full,lapserate_corr,entrainment_profile)
end

# use Grid in Parameters if not provided
Geometry(P::Parameters) = Geometry(P,P.Grid)

"""
    S = SpectralTransform(P::Parameters)

Generator function for a SpectralTransform struct pulling in parameters from a Parameters struct."""
function SpeedyTransforms.SpectralTransform(P::Parameters)
    (;NF, Grid, trunc, dealiasing, recompute_legendre, legendre_shortcut) = P
    return SpectralTransform(NF,Grid,trunc,recompute_legendre;legendre_shortcut,dealiasing)
end