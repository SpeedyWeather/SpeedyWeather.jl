"""
    HD = HorizontalDiffusion(...)

Horizontal Diffusion struct containing all the preallocated arrays for the calculation of horizontal diffusion
and orographic correction for temperature and humidity.

    # Explicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping::Matrix{NF}                 # for temperature and vorticity (explicit)
    damping_div::Matrix{NF}             # for divergence (explicit)
    damping_strat::Matrix{NF}           # for extra diffusion in the stratosphere (explicit)
    
    # Implicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping_impl::Matrix{NF}            # for temperature and vorticity (implicit)
    damping_div_impl::Matrix{NF}        # for divergence (implicit)
    damping_strat_impl::Matrix{NF}      # for extra diffusion in the stratosphere (implicit)
    
    # Vertical component of orographic correction
    temp_correction_vert::Vector{NF}    # for temperature
    humid_correction_vert::Vector{NF}   # for humidity
    
    # Horizontal component of orographic correction (in spectral space)
    temp_correction_horizontal::Matrix{Complex{NF}}     # for temperature
    humid_correction_horizontal::Matrix{Complex{NF}}    # for humidity.
"""
struct HorizontalDiffusion{NF<:AbstractFloat}   # Number format NF
    # Explicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping::Matrix{NF}                 # for temperature and vorticity (explicit)
    damping_div::Matrix{NF}             # for divergence (explicit)
    damping_strat::Matrix{NF}           # for extra diffusion in the stratosphere (explicit)
    
    # Implicit part of the diffusion, precalculated damping coefficients for each spectral mode
    damping_impl::Matrix{NF}            # for temperature and vorticity (implicit)
    damping_div_impl::Matrix{NF}        # for divergence (implicit)
    damping_strat_impl::Matrix{NF}      # for extra diffusion in the stratosphere (implicit)
    
    # Vertical component of orographic correction
    temp_correction_vert::Vector{NF}    # for temperature
    humid_correction_vert::Vector{NF}   # for humidity
    
    # Horizontal component of orographic correction (in spectral space)
    temp_correction_horizontal::Matrix{Complex{NF}}     # for temperature
    humid_correction_horizontal::Matrix{Complex{NF}}    # for humidity
end

"""
    HD = HorizontalDiffusion(::Parameters,::GeoSpectral,::Boundaries)

Generator function for a HorizontalDiffusion struct `HD`. Precalculates damping matrices for
horizontal hyperdiffusion for temperature, vorticity and divergence, with an implicit term
and an explicit term. Also precalculates correction terms (horizontal and vertical) for
temperature and humidity.
"""
function HorizontalDiffusion(   P::Parameters,      # Parameter struct
                                C::Constants,       # Constants struct
                                G::GeoSpectral,     # Geometry and spectral struct
                                B::Boundaries)      # Boundaries struct

    # DIFFUSION
    @unpack lmax,mmax,radius = G.spectral_transform
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
    damping = zeros(lmax+1,mmax+1)              # for temperature and vorticity (explicit)
    damping_div = zeros(lmax+1,mmax+1)          # for divergence (explicit)
    damping_strat = zeros(lmax+1,mmax+1)        # for extra diffusion in the stratosphere (explicit)

    # Damping coefficients for implicit part of the diffusion (= 1/(1+2Δtν∇²ⁿ))
    damping_impl = zeros(lmax+1,mmax+1)         # for temperature and vorticity (implicit)
    damping_div_impl = zeros(lmax+1,mmax+1)     # for divergence (implicit)
    damping_strat_impl = zeros(lmax+1,mmax+1)   # for extra diffusion in the stratosphere (implicit)

    # PRECALCULATE the damping coefficients for every spectral mode
    for m in 1:mmax+1                           # fill only the lower triangle
        for l in m:lmax+1
            # eigenvalue is l*(l+1), but 1-based here l→l-1
            norm_eigenvalue = l*(l-1)/largest_eigenvalue        # normal diffusion ∇²
            norm_eigenvalueⁿ = norm_eigenvalue^diffusion_power  # hyper diffusion ∇²ⁿ

            # Explicit part (=ν∇²ⁿ)
            # convert diffusion time scales to damping frequencies [1/s] times norm. eigenvalue
            damping[l,m] = norm_eigenvalueⁿ/(3600*diffusion_time)*radius                # temperature/vorticity
            damping_div[l,m] = norm_eigenvalueⁿ/(3600*diffusion_time_div)*radius        # divergence
            damping_strat[l,m] = norm_eigenvalue/(3600*diffusion_time_strat)*radius     # stratosphere (no hyperdiff)

            # and implicit part of the diffusion (= 1/(1+2Δtν∇²ⁿ))
            damping_impl[l,m] = 1/(1+2Δt*damping[l,m])                 # for temperature/vorticity
            damping_div_impl[l,m] = 1/(1+2Δt*damping_div[l,m])         # for divergence
            damping_strat_impl[l,m] = 1/(1+2Δt*damping_strat[l,m])     # for stratosphere (no hyperdiffusion)
        end
    end

    # OROGRAPHIC CORRECTION
    @unpack nlon, nlat, nlev, σ_levels_full = G.geometry
    @unpack gravity, R_gas, lapse_rate, scale_height, scale_height_humid, relhumid_ref = P
    @unpack geopot_surf_grid = B    #TODO is currently not contained in the Boundaries struct B

    # Orographic correction terms for temperature and humidity (vertical component)
    lapse_rate_gravity = lapse_rate/(1000gravity)       # lapse rate in [K/km] convert to [K/m] with /1000
    R_lapse_rate_gravity = R_gas*lapse_rate_gravity     
    scale_height_ratio = scale_height/scale_height_humid

    # preallocate (high precision, conversion to NF later)
    temp_correction_vert = zeros(nlev)      # Vertical component of orographic correction for temperature
    humid_correction_vert = zeros(nlev)     # Vertical component of orographic correction for humidity

    for k in 2:nlev
        temp_correction_vert[k] = σ_levels_full[k]^R_lapse_rate_gravity
        if k > 2
            humid_correction_vert[k] = σ_levels_full[k]^scale_height_ratio
        end
    end

    # Orographic correction term for temperature (horizontal component)
    horizontal_correction = zeros(nlon, nlat)       # in grid-point space

    for j in 1:nlat
        for i = 1:nlon
            horizontal_correction[i,j] = lapse_rate_gravity*geopot_surf_grid[i,j]
        end
    end

    # transform correction to spectral space
    temp_correction_horizontal = spectral(horizontal_correction,one_more_l=true)

    # Orographic correction terms for humidity (horizontal component)
    horizontal_correction .= relhumid_ref           # relative humidity reference value
    humid_correction_horizontal = spectral(horizontal_correction,one_more_l=true)

    # convert to number format NF here
    return HorizontalDiffusion{P.NF}(   damping,damping_div,damping_strat,
                                        damping_impl,damping_div_impl,damping_strat_impl,
                                        temp_correction_vert,humid_correction_vert,
                                        temp_correction_horizontal,humid_correction_horizontal)
end

"""
    horizontal_diffusion!(  tendency::AbstractMatrix{Complex{NF}}, # tendency of a 
                            A::AbstractMatrix{Complex{NF}},        # spectral horizontal field
                            damp_expl::AbstractMatrix{NF},         # explicit spectral damping
                            damp_impl::AbstractMatrix{NF}          # implicit spectral damping
                            ) where {NF<:AbstractFloat}

Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `damp_expl` and an implicit part `damp_impl`, such that both can be calculated
in a single forward step by using `A` as well as its tendency `tendency`."""
function horizontal_diffusion!( tendency::AbstractMatrix{Complex{NF}}, # tendency of a 
                                A::AbstractMatrix{Complex{NF}},        # spectral horizontal field
                                damp_expl::AbstractMatrix{NF},         # explicit spectral damping
                                damp_impl::AbstractMatrix{NF}          # implicit spectral damping
                                ) where {NF<:AbstractFloat}

    lmax,mmax = size(A) .- 1            # degree l, order m but 0-based
    @boundscheck size(A) == size(tendency) || throw(BoundsError())
    @boundscheck size(A) == size(damp_expl) || throw(BoundsError())
    @boundscheck size(A) == size(damp_impl) || throw(BoundsError())
    
    @inbounds for m in 1:mmax+1         # loop through all spectral modes 
        for l in m:lmax+1
            tendency[l,m] = (tendency[l,m] - damp_expl[l,m]*A[l,m])*damp_impl[l,m]
        end
    end
end

"""
    horizontal_diffusion!(  tendency::AbstractArray{Complex{NF},3}, # tendency of a
                            A::AbstractArray{Complex{NF},3},        # spectral spatial field
                            damp_expl::AbstractMatrix{NF},          # explicit spectral damping
                            damp_impl::AbstractMatrix{NF}           # implicit spectral damping
                            ) where {NF<:AbstractFloat}             # number format NF

Apply horizontal diffusion to a 3D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The diffusion is applied layer by layer. The implicit
diffusion of the next time step is split into an explicit part `damp_expl` and an implicit part
`damp_impl`, such that both can be calculated in a single forward step by using `A` as well as
its tendency `tendency`."""
function horizontal_diffusion!( tendency::AbstractArray{Complex{NF},3}, # tendency of a
                                A::AbstractArray{Complex{NF},3},        # spectral spatial field
                                damp_expl::AbstractMatrix{NF},          # explicit spectral damping
                                damp_impl::AbstractMatrix{NF}           # implicit spectral damping
                                ) where {NF<:AbstractFloat}             # number format NF
    _,_,nlev = size(A)
    @boundscheck size(A) == size(tendency) || throw(BoundsError())
    
    for k in 1:nlev
        A_layer = view(A,:,:,k)
        tendency_layer = view(tendency,:,:,k)
        horizontal_diffusion!(tendency_layer, A_layer, damp_expl, damp_impl)
    end
end

"""
    stratospheric_zonal_drag!(  tendency::AbstractArray{Complex{NF},3}, # tendency of
                                A::AbstractArray{Complex{NF},3},        # spectral vorticity or divergence
                                drag::Real                              # drag coefficient [1/s]
                                ) where {NF<:AbstractFloat}             # number format NF

Zonal drag in the uppermost layer of the stratosphere of 3D spectral field `A` (vorticity or divergence).
Drag is applied explicitly to the time step in `A` and its tendency `tendency` is changed in-place.
`drag` is the drag coefficient of unit 1/s.
"""
function stratospheric_zonal_drag!( tendency::AbstractArray{Complex{NF},3}, # tendency of
                                    A::AbstractArray{Complex{NF},3},        # spectral vorticity or divergence
                                    drag::Real                              # drag coefficient [1/s]
                                    ) where {NF<:AbstractFloat}             # number format NF
    
    lmax,mmax,nlev = size(A)    # spherical harmonic degree l, order m, number of vertical levels nlev
    lmax -= 1                   # convert to 0-based l,m
    mmax -= 1
    @boundscheck size(A) == size(tendency) || throw(BoundsError())

    drag_NF = convert(NF,drag)

    @inbounds for l in 1:lmax+1     # loop over degree l, but 1-based
        # size(A) = lmax x mmax x nlev, nlev = 1 is uppermost model level
        # apply drag only to largest zonal wavenumber (m = 0) and in the uppermost model level (k=1)
        tendency[l,1,1] = tendency[l,1,1] - drag_NF*A[l,1,1]
    end
end

"""Orographic temperature correction for absolute temperature to be applied before the horizontal diffusion."""
function orographic_correction!(A_corrected::AbstractArray{Complex{NF},3},  # correction of 
                                A::AbstractArray{Complex{NF},3},            # 3D spectral temperature or humidity
                                correction_horizontal::Matrix{Complex{NF}}, # horizontal correction matrix
                                correction_vertical::Vector{NF},            # vertical correction vector
                                ) where NF
    
    lmax,mmax,nlev = size(A)    # degree l, order m of the spherical harmonics
    lmax -= 1                   # convert to 0-based
    mmax -= 1

    @boundscheck size(A) == size(A_corrected) || throw(BoundsError())
    @boundscheck (lmax+1,mmax+1) == size(correction_horizontal) || throw(BoundsError())
    @boundscheck (nlev,) == size(correction_vertical) || throw(BoundsError())

    @inbounds for k in 1:nlev       # vertical levels
        for m in 1:mmax+1           # order of spherical harmonics
            for l in m:lmax+1       # degree of spherical harmonics
                A_corrected[l,m,k] = A[l,m,k] + hori_correction[l,m]*vert_correction[k]
            end
        end
    end
end