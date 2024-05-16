abstract type AbstractHorizontalDiffusion <: AbstractModelComponent end

export HyperDiffusion

"""
Horizontal hyper diffusion of vor, div, temp, humid; implicitly in spectral space
with a `power` of the Laplacian (default = 4) and the strength controlled by
`time_scale` (default = 1 hour). By default, the time scale (=1/strength of diffusion)
is reduced with increasing resolution, and the power is linearly decreased in the vertical
above the `tapering_σ` sigma level to `power_stratosphere` (default 2). Fields and options are
$(TYPEDFIELDS)"""
@kwdef mutable struct HyperDiffusion{NF} <: AbstractHorizontalDiffusion
    # DIMENSIONS
    "spectral resolution"
    trunc::Int

    "number of vertical levels"
    nlev::Int
    
    # PARAMETERS
    "[OPTION] power of Laplacian"
    power::Float64 = 4.0
    
    "[OPTION] diffusion time scale"
    time_scale::Second = Minute(60)
    
    "[OPTION] stronger diffusion with resolution? 0: constant with trunc, 1: (inverse) linear with trunc, etc"
    resolution_scaling::Float64 = 0.5

    # incrased diffusion in stratosphere
    "[OPTION] different power for tropopause/stratosphere"
    power_stratosphere::Float64 = 2.0
    
    "[OPTION] linearly scale towards power_stratosphere above this σ"
    tapering_σ::Float64 = 0.2

    # ARRAYS, precalculated for each spherical harmonics degree and vertical layer
    ∇²ⁿ::Vector{Vector{NF}} = [zeros(NF, trunc+2) for _ in 1:nlev]           # explicit part
    ∇²ⁿ_implicit::Vector{Vector{NF}} = [ones(NF, trunc+2) for _ in 1:nlev]   # implicit part
end

"""$(TYPEDSIGNATURES)
Generator function based on the resolutin in `spectral_grid`.
Passes on keyword arguments."""
function HyperDiffusion(spectral_grid::SpectralGrid; kwargs...)
    (; NF, trunc, nlev) = spectral_grid        # take resolution parameters from spectral_grid
    return HyperDiffusion{NF}(; trunc, nlev, kwargs...)
end

"""$(TYPEDSIGNATURES)
Precomputes the hyper diffusion terms in `scheme` based on the
model time step, and possibly with a changing strength/power in
the vertical.
"""
function initialize!(   scheme::HyperDiffusion,
                        model::ModelSetup)
    initialize!(scheme, model.geometry, model.time_stepping)
end

"""$(TYPEDSIGNATURES)
Precomputes the hyper diffusion terms for all layers based on the
model time step in `L`, the vertical level sigma level in `G`."""
function initialize!(   
    scheme::HyperDiffusion,
    G::AbstractGeometry,
    L::AbstractTimeStepper,
)
    (; trunc, nlev, resolution_scaling, ∇²ⁿ, ∇²ⁿ_implicit) = scheme
    (; power, power_stratosphere, tapering_σ) = scheme
    (; Δt, radius) = L

    # Reduce diffusion time scale (=increase diffusion, always in seconds) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    time_scale = scheme.time_scale.value/radius * (32/(trunc+1))^resolution_scaling

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -trunc*(trunc+1)
    
    for k in 1:nlev
        # VERTICAL TAPERING for the stratosphere
        # go from 1 to 0 between σ=0 and tapering_σ
        σ = G.σ_levels_full[k]
        tapering = max(0, (tapering_σ-σ)/tapering_σ)         # ∈ [0, 1]
        p = power + tapering*(power_stratosphere - power)

        for l in 0:trunc    # diffusion for every degree l, but indendent of order m
            eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

            # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
            ∇²ⁿ[k][l+1] = -eigenvalue_norm^p/time_scale
            
            # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
            ∇²ⁿ_implicit[k][l+1] = 1/(1-2Δt*∇²ⁿ[k][l+1])           
        end
        
        # last degree is only used by vector quantities; set to zero for implicit and explicit
        # to set any tendency at lmax+1,1:mmax to zero (what it should be anyway)
        ∇²ⁿ[k][trunc+2] = 0
        ∇²ⁿ_implicit[k][trunc+2] = 0
    end
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to a 2D field `A` in spectral space by updating its tendency `tendency`
with an implicitly calculated diffusion term. The implicit diffusion of the next time step is split
into an explicit part `∇²ⁿ_expl` and an implicit part `∇²ⁿ_impl`, such that both can be calculated
in a single forward step by using `A` as well as its tendency `tendency`."""
function horizontal_diffusion!( tendency::LowerTriangularMatrix{Complex{NF}},   # tendency of a 
                                A::LowerTriangularMatrix{Complex{NF}},          # spectral horizontal field
                                ∇²ⁿ_expl::AbstractVector{NF},                   # explicit spectral damping
                                ∇²ⁿ_impl::AbstractVector{NF}                    # implicit spectral damping
                                ) where NF
    lmax, mmax = size(tendency)      # 1-based
    @boundscheck size(tendency) == size(A) || throw(BoundsError)
    @boundscheck lmax <= length(∇²ⁿ_expl) == length(∇²ⁿ_impl) || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax   # loops over all columns/order m
        for l in m:lmax-1       # but skips the lmax+2 degree (1-based)
            lm += 1             # single index lm corresponding to harmonic l, m
            tendency[lm] = (tendency[lm] + ∇²ⁿ_expl[l]*A[lm])*∇²ⁿ_impl[l]
        end
        lm += 1                 # skip last row for scalar quantities
    end
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity in the BarotropicModel."""
function horizontal_diffusion!( diagn::DiagnosticVariablesLayer,
                                progn::PrognosticLayerTimesteps,
                                model::Barotropic,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    ∇²ⁿ = model.horizontal_diffusion.∇²ⁿ[diagn.k]
    ∇²ⁿ_implicit = model.horizontal_diffusion.∇²ⁿ_implicit[diagn.k]

    # Barotropic model diffuses vorticity (only variable)
    (; vor) = progn.timesteps[lf]
    (; vor_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity and divergence in the ShallowWaterModel."""
function horizontal_diffusion!( progn::PrognosticLayerTimesteps,
                                diagn::DiagnosticVariablesLayer,
                                model::ShallowWater,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    ∇²ⁿ = model.horizontal_diffusion.∇²ⁿ[diagn.k]
    ∇²ⁿ_implicit = model.horizontal_diffusion.∇²ⁿ_implicit[diagn.k]

    # ShallowWater model diffuses vorticity and divergence
    (; vor, div) = progn.timesteps[lf]
    (; vor_tend, div_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(div_tend, div, ∇²ⁿ, ∇²ⁿ_implicit)
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion applied to vorticity, divergence, temperature, and
humidity (PrimitiveWet only) in the PrimitiveEquation models."""
function horizontal_diffusion!( progn::PrognosticLayerTimesteps,
                                diagn::DiagnosticVariablesLayer,
                                model::PrimitiveEquation,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    ∇²ⁿ = model.horizontal_diffusion.∇²ⁿ[diagn.k]
    ∇²ⁿ_implicit = model.horizontal_diffusion.∇²ⁿ_implicit[diagn.k]

    # Primitive equation models diffuse vor, divergence, temp (and humidity for wet core)
    (; vor, div, temp, humid) = progn.timesteps[lf]
    (; vor_tend, div_tend, temp_tend, humid_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(div_tend, div, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(temp_tend, temp, ∇²ⁿ, ∇²ⁿ_implicit)
    model isa PrimitiveWet && horizontal_diffusion!(humid_tend, humid, ∇²ⁿ, ∇²ⁿ_implicit)
end