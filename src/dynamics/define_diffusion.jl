Base.@kwdef struct HyperDiffusion <: DiffusionParameters
    # hyperdiffusion for temp, vor, div everywhere
    # several powers of Laplacians, default 4 and 2, are added
    # with respective time scales and scalings with resolution
    power::Float64 = 4.0                # Powers of Laplacians
    time_scale::Float64 = 2.4           # Diffusion time scales [hrs]
    resolution_scaling::Float64 = 1.0   # 1: (inverse) linear with T
                                        # 2: (inverse) quadratic, etc

    # additional diffusion in stratosphere
    power_stratosphere::Float64 = 2.0   # different power for stratosphere
    tapering_σ::Float64 = 0.2           # scale towards that power linearly above this σ

    # increase diffusion based on high vorticity levels
    adaptive::Bool = true               # swith on/off
    vor_max::Float64 = 1e-4             # [1/s] above this, diffusion is increased
    adaptive_strength::Float64 = 2.0    # increase strength above vor_max by this factor
                                        # times max(abs(vor))/vor_max
end

""" 
    HD = HorizontalDiffusion(...)

Horizontal Diffusion struct containing all the preallocated arrays for the calculation
of horizontal diffusion."""
struct HorizontalDiffusion{NF<:AbstractFloat}   # Number format NF
    
    lmax::Int               # max degree l of spherical harmonics
    time_scale::NF          # strength of diffusion ~1/time_scale

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
    @unpack power, time_scale, resolution_scaling = scheme
    @unpack power_stratosphere, tapering_σ = scheme
    @unpack Δt = C
    @unpack nlev, σ_levels_full = G

    # Reduce diffusion time scale (=increase diffusion) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = 1/radius*3600*scheme.time_scale * (32/(lmax+1))^resolution_scaling

    # PREALLOCATE as vector as only dependend on degree l
    # Damping coefficients for explicit part of the diffusion (=ν∇²ⁿ)
    ∇²ⁿ = [zeros(NF,lmax+2) for _ in 1:nlev]                # for temp, vor, div (explicit)
    ∇²ⁿ_implicit = [zeros(NF,lmax+2) for _ in 1:nlev]       # Implicit part (= 1/(1+2Δtν∇²ⁿ))
    HD = HorizontalDiffusion{NF}(lmax,time_scale,∇²ⁿ,∇²ⁿ_implicit)  # initialize struct

    # PRECOMPUTE diffusion operators based on non-adaptive diffusion
    vor_max = 0                                 # not increased diffusion initially
    adapt_diffusion!(HD,scheme,vor_max,G,C)     # fill diffusion arrays
    return HD
end

function adapt_diffusion!(  HD::HorizontalDiffusion,
                            scheme::HyperDiffusion,
                            vor_max::Real,
                            G::Geometry,
                            C::DynamicsConstants)
    nlev = length(HD.∇²ⁿ)
    for k in 1:nlev
        adapt_diffusion!(HD,scheme,vor_max,k,G,C)
    end
end

function adapt_diffusion!(  HD::HorizontalDiffusion,
                            scheme::HyperDiffusion,
                            vor_max::Real,
                            k::Int,
                            G::Geometry,
                            C::DynamicsConstants)

    (;lmax,time_scale,∇²ⁿ,∇²ⁿ_implicit) = HD
    (;power, power_stratosphere, tapering_σ) = scheme
    (;Δt) = C
    σ = G.σ_levels_full[k]

    # ADAPTIVE/FLOW AWARE
    # increase diffusion if maximum vorticity per layer is larger than scheme.vor_max
    if scheme.adaptive
        time_scale /= max(1,scheme.adaptive_strength*vor_max/scheme.vor_max)
    end

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -lmax*(lmax+1)
    
    # VERTICAL TAPERING for the stratosphere
    # go from 1 to 0 between σ=0 and tapering_σ
    tapering = max(0,(tapering_σ-σ)/tapering_σ)         # ∈ [0,1]
    p = power + tapering*(power_stratosphere - power) 
        
    @inbounds for l in 0:lmax+1   # diffusion for every degree l, but indendent of order m
        eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

        # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
        ∇²ⁿ[k][l+1] = -eigenvalue_norm^p/time_scale
        
        # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
        ∇²ⁿ_implicit[k][l+1] = 1/(1-2Δt*∇²ⁿ[k][l+1])           
    end
end