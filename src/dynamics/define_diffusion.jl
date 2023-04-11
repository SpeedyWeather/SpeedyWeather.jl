Base.@kwdef struct HyperDiffusion <: DiffusionParameters
    # hyperdiffusion for temp, vor, div everywhere
    power::Int = 4                              # Power n of Laplacian in horizontal diffusion ∇²ⁿ
    time_scale::Float64 = 2.5                   # Diffusion time scale [hrs] for temp, vor, div
    
    # additional diffusion in stratosphere
    power_stratosphere::Int = 2                 # different power for stratosphere
    tapering_σ::Float64 = 0.2                   # scale towards that power linearly above this σ

    # reduce time scale of diffusion linearly with increasing model resolution?
    scale_with_resolution::Bool = true          # e.g. 2.5hrs for T31, 1.25 hrs for T63 etc
end

""" 
    HD = HorizontalDiffusion(...)

Horizontal Diffusion struct containing all the preallocated arrays for the calculation
of horizontal diffusion."""
struct HorizontalDiffusion{NF<:AbstractFloat}   # Number format NF
    
    # (Hyper) diffusion, precalculated for each spherical harm degree
    # and for each layer (to allow for varying orders/strength in the vertical)
    ∇²ⁿ::Vector{Vector{NF}}                     # explicit part
    ∇²ⁿ_implicit::Vector{Vector{NF}}            # implicit part
end

"""
    HD = HorizontalDiffusion(::Parameters,::GeoSpectral,::Boundaries)

Generator function for a HorizontalDiffusion struct `HD`. Precalculates damping matrices for
horizontal hyperdiffusion for temperature, vorticity and divergence, with an implicit term
and an explicit term. Also precalculates correction terms (horizontal and vertical) for
temperature and humidity.
"""
function HorizontalDiffusion(   scheme::HyperDiffusion,
                                P::Parameters,
                                C::DynamicsConstants,
                                G::Geometry,
                                S::SpectralTransform{NF}) where NF
    @unpack lmax,mmax = S
    @unpack radius = P.planet
    @unpack power, power_stratosphere, tapering_σ = scheme
    @unpack Δt = C
    @unpack nlev, σ_levels_full = G

    # reuce diffusion time scale (=increase diffusion) with resolution?
    if scheme.scale_with_resolution
        # use values in scheme for T31 (=32 here) and decrease linearly with lmax+1
        time_scale = scheme.time_scale * (32/(lmax+1))
    else
        # @unpack time_scale, time_scale_stratosphere = scheme
        @unpack time_scale = scheme
    end

    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -lmax*(lmax+1)

    # PREALLOCATE as vector as only dependend on degree l
    # Damping coefficients for explicit part of the diffusion (=ν∇²ⁿ)
    ∇²ⁿ = [zeros(NF,lmax+2) for _ in 1:nlev]                # for temp, vor, div (explicit)
    ∇²ⁿ_implicit = [zeros(NF,lmax+2) for _ in 1:nlev]       # Implicit part (= 1/(1+2Δtν∇²ⁿ))

    for k in 1:nlev         # every layer is a combination of two Laplacians of different orders

        # tapering: go from 1 to 0 between σ=0 and tapering_σ
        σ = σ_levels_full[k]
        tapering = max(0,(tapering_σ-σ)/tapering_σ)         # ∈ [0,1]
        p = power + tapering*(power_stratosphere - power)     
        
        for l in 0:lmax+1   # PRECALCULATE for every degree l, but indendent of order m
            eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

            # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
            # time scale [hrs] *3600-> [s] and add additional diffusion for stratosphere
            ∇²ⁿ[k][l+1] = -eigenvalue_norm^p/(3600*time_scale)*radius
            
            # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
            ∇²ⁿ_implicit[k][l+1] = 1/(1-2Δt*∇²ⁿ[k][l+1])           
        end
    end

    return HorizontalDiffusion(∇²ⁿ,∇²ⁿ_implicit)
end