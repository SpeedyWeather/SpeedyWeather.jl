Base.@kwdef struct HyperDiffusion <: DiffusionParameters
    # hyperdiffusion for temp, vor, div everywhere
    # several powers of Laplacians, default 4 and 2, are added
    # with respective time scales and scalings with resolution
    powers::Vector{Float64} = [4.0,2.0]                 # Powers of Laplacians
    time_scales::Vector{Float64} = [2.4,12.0]           # Diffusion time scales [hrs]
    resolution_scalings::Vector{Float64} = [1.0,2.0]    # 1: (inverse) linear with T
                                                        # 2: (inverse) quadratic, etc
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
    @unpack powers, resolution_scalings = scheme
    @unpack Δt = C
    @unpack nlev = G

    # Reduce diffusion time scale (=increase diffusion) with resolution
    time_scales = zero(scheme.time_scales)
    for i in eachindex(time_scales,resolution_scalings)
        # use values in scheme for T31 (=32 here) and decrease with lmax+1
        # time scale [hrs] *3600-> [s]
        time_scales[i] = 3600*scheme.time_scales[i] * (32/(lmax+1))^resolution_scalings[i]
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

    for k in 1:nlev         # every layer is a combination of Laplacians of different orders
        for l in 0:lmax+1   # PRECALCULATE for every degree l, but indendent of order m
            eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

            # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
            ∇²ⁿ[k][l+1] = 0
            for i in eachindex(powers,time_scales)
                power = powers[i]
                time_scale = time_scales[i]
                ∇²ⁿ[k][l+1] += -eigenvalue_norm^power/time_scale*radius
            end
            
            # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
            ∇²ⁿ_implicit[k][l+1] = 1/(1-2Δt*∇²ⁿ[k][l+1])           
        end
    end

    return HorizontalDiffusion(∇²ⁿ,∇²ⁿ_implicit)
end