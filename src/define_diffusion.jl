"""
    HD = HorizontalDiffusion(...)

Horizontal Diffusion struct containing all the preallocated arrays for the calculation of horizontal diffusion
and orographic correction for temperature and humidity.
"""
struct HorizontalDiffusion{NF<:AbstractFloat}       # Number format NF
    # Explicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping::LowerTriangularMatrix{NF}              # for temperature and vorticity (explicit)
    damping_div::LowerTriangularMatrix{NF}          # for divergence (explicit)
    damping_strat::LowerTriangularMatrix{NF}        # for extra diffusion in the stratosphere (explicit)
    
    # Implicit part of LowerTriangular diffusion, precalculated damping coefficients for each spectral mode
    damping_impl::LowerTriangularMatrix{NF}         # for temperature and vorticity (implicit)
    damping_div_impl::LowerTriangularMatrix{NF}     # for divergence (implicit)
    damping_strat_impl::LowerTriangularMatrix{NF}   # for extra diffusion in the stratosphere (implicit)
    
    # Vertical component of orographic correction
    temp_correction_vert::Vector{NF}                # for temperature
    humid_correction_vert::Vector{NF}               # for humidity
    
    # Horizontal component of orographic correction (in spectral space)
    temp_correction_horizontal::LowerTriangularMatrix{Complex{NF}}
    humid_correction_horizontal::LowerTriangularMatrix{Complex{NF}}
end

"""
    HD = HorizontalDiffusion(::Parameters,::GeoSpectral,::Boundaries)

Generator function for a HorizontalDiffusion struct `HD`. Precalculates damping matrices for
horizontal hyperdiffusion for temperature, vorticity and divergence, with an implicit term
and an explicit term. Also precalculates correction terms (horizontal and vertical) for
temperature and humidity.
"""
function HorizontalDiffusion(   P::Parameters,          # Parameter struct
                                C::Constants,           # Constants struct
                                G::Geometry,            # Geometry struct
                                S::SpectralTransform,   # SpectralTransform struct 
                                B::Boundaries)          # Boundaries struct

    # DIFFUSION
    @unpack lmax,mmax = S
    @unpack radius_earth = G
    @unpack diffusion_power, diffusion_time, diffusion_time_div = P
    @unpack diffusion_time_strat, damping_time_strat = P
    @unpack Δt = C

    # Diffusion is applied by multiplication of the (absolute) eigenvalues of the Laplacian l*(l+1)
    # normalise by the largest eigenvalue lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the time scale diffusion_time
    # raise to a power of the Laplacian for hyperdiffusion (=less damping for smaller wavenumbers)
    largest_eigenvalue = lmax*(lmax+1)

    # PREALLOCATE
    # conversion to number format NF later, one more degree l for meridional gradient recursion
    # Damping coefficients for explicit part of the diffusion (=ν∇²ⁿ)
    # while precalculated for spectral space, store only the real part as entries are real
    LTM = LowerTriangularMatrix
    damping = zeros(LTM,lmax+2,mmax+1)              # for temperature and vorticity (explicit)
    damping_div = zeros(LTM,lmax+2,mmax+1)          # for divergence (explicit)
    damping_strat = zeros(LTM,lmax+2,mmax+1)        # for extra diffusion in the stratosphere (explicit)

    # Damping coefficients for implicit part of the diffusion (= 1/(1+2Δtν∇²ⁿ))
    damping_impl = zeros(LTM,lmax+2,mmax+1)         # for temperature and vorticity (implicit)
    damping_div_impl = zeros(LTM,lmax+2,mmax+1)     # for divergence (implicit)
    damping_strat_impl = zeros(LTM,lmax+2,mmax+1)   # for extra diffusion in the stratosphere (implicit)

    # PRECALCULATE the damping coefficients for every spectral mode
    R = radius_earth                                # convenience
    for m in 1:mmax+1                               # fill only the lower triangle
        for l in m:lmax+1
            # eigenvalue is l*(l+1), but 1-based here l→l-1
            norm_eigenvalue = l*(l-1)/largest_eigenvalue        # normal diffusion ∇²
            norm_eigenvalueⁿ = norm_eigenvalue^diffusion_power  # hyper diffusion ∇²ⁿ

            # Explicit part (=ν∇²ⁿ)
            # convert diffusion time scales to damping frequencies [1/s] times norm. eigenvalue
            damping[l,m] = norm_eigenvalueⁿ/(3600*diffusion_time)*R                # temperature/vorticity
            damping_div[l,m] = norm_eigenvalueⁿ/(3600*diffusion_time_div)*R        # divergence
            damping_strat[l,m] = norm_eigenvalue/(3600*diffusion_time_strat)*R     # stratosphere (no hyperdiff)

            # and implicit part of the diffusion (= 1/(1+2Δtν∇²ⁿ))
            damping_impl[l,m] = 1/(1+2Δt*damping[l,m])                 # for temperature/vorticity
            damping_div_impl[l,m] = 1/(1+2Δt*damping_div[l,m])         # for divergence
            damping_strat_impl[l,m] = 1/(1+2Δt*damping_strat[l,m])     # for stratosphere (no hyperdiffusion)
        end
    end

    if P.model <: Barotropic || P.model <: ShallowWater                 # orographic correction not needed
        
        temp_correction_vert        = zeros(0)                          # create dummy arrays
        humid_correction_vert       = zeros(0)
        temp_correction_horizontal  = zeros(LowerTriangularMatrix{Complex{P.NF}},0,0) 
        humid_correction_horizontal = zeros(LowerTriangularMatrix{Complex{P.NF}},0,0) 

    else    # P.model <: PrimitiveEquation, orographic correction only needed then

        temp_correction_vert        = zeros(0)                          # create dummy arrays
        humid_correction_vert       = zeros(0)
        temp_correction_horizontal  = zeros(LowerTriangularMatrix{Complex{P.NF}},0,0) 
        humid_correction_horizontal = zeros(LowerTriangularMatrix{Complex{P.NF}},0,0) 

        # # OROGRAPHIC CORRECTION
        # @unpack nlon, nlat, nlev, σ_levels_full = G
        # @unpack gravity, R_gas, lapse_rate, scale_height, scale_height_humid, relhumid_ref = P
        # @unpack geopot_surf_grid = B    #TODO is currently not contained in the Boundaries struct B

        # # Orographic correction terms for temperature and humidity (vertical component)
        # lapse_rate_gravity = lapse_rate/(1000gravity)       # lapse rate in [K/km] convert to [K/m] with /1000
        # R_lapse_rate_gravity = R_gas*lapse_rate_gravity     
        # scale_height_ratio = scale_height/scale_height_humid

        # # preallocate (high precision, conversion to NF later)
        # temp_correction_vert = zeros(nlev)      # Vertical component of orographic correction for temperature
        # humid_correction_vert = zeros(nlev)     # Vertical component of orographic correction for humidity

        # for k in 2:nlev
        #     temp_correction_vert[k] = σ_levels_full[k]^R_lapse_rate_gravity
        #     if k > 2
        #         humid_correction_vert[k] = σ_levels_full[k]^scale_height_ratio
        #     end
        # end

        # # Orographic correction term for temperature (horizontal component)
        # horizontal_correction = zeros(nlon, nlat)       # in grid-point space

        # for j in 1:nlat
        #     for i = 1:nlon
        #         horizontal_correction[i,j] = lapse_rate_gravity*geopot_surf_grid[i,j]
        #     end
        # end

        # # transform correction to spectral space
        # temp_correction_horizontal = spectral(horizontal_correction,one_more_l=true)

        # # Orographic correction terms for humidity (horizontal component)
        # horizontal_correction .= relhumid_ref           # relative humidity reference value
        # humid_correction_horizontal = spectral(horizontal_correction,one_more_l=true)
    end

    # convert to number format NF here
    return HorizontalDiffusion{P.NF}(   damping,damping_div,damping_strat,
                                        damping_impl,damping_div_impl,damping_strat_impl,
                                        temp_correction_vert,humid_correction_vert,
                                        temp_correction_horizontal,humid_correction_horizontal)
end