abstract type AbstractHorizontalDiffusion <: AbstractModelComponent end

export HyperDiffusion

"""
Struct for horizontal hyper diffusion of vor, div, temp; implicitly in spectral space
with a `power` of the Laplacian (default=4) and the strength controlled by
`time_scale`. Options exist to scale the diffusion by resolution, and adaptive
depending on the current vorticity maximum to increase diffusion in active
layers. Furthermore the power can be decreased above the `tapering_σ` to
`power_stratosphere` (default 2). For Barotropic, ShallowWater,
the default non-adaptive constant-time scale hyper diffusion is used. Options are
$(TYPEDFIELDS)"""
Base.@kwdef mutable struct HyperDiffusion{NF} <: AbstractHorizontalDiffusion
    # DIMENSIONS
    "spectral resolution"
    trunc::Int

    "number of vertical levels"
    nlev::Int
    
    # PARAMETERS
    "power of Laplacian"
    power::Float64 = 4.0
    
    "diffusion time scale"
    time_scale::Second = Minute(144)
    
    "stronger diffusion with resolution? 0: constant with trunc, 1: (inverse) linear with trunc, etc"
    resolution_scaling::Float64 = 0.5

    # incrased diffusion in stratosphere
    "different power for tropopause/stratosphere"
    power_stratosphere::Float64 = 2.0
    
    "linearly scale towards power_stratosphere above this σ"
    tapering_σ::Float64 = 0.2

    # increase diffusion based on high vorticity levels
    "adaptive = higher diffusion for layers with higher vorticity levels."
    adaptive::Bool = true               # swith on/off
    
    "above this (absolute) vorticity level [1/s], diffusion is increased"
    vor_max::Float64 = 1e-4
    
    "increase strength above `vor_max` by this factor times `max(abs(vor))/vor_max`"
    adaptive_strength::Float64 = 2.0

    # ARRAYS, precalculated for each spherical harmonics degree
    # Barotropic and ShallowWater are fine with a constant time scale 
    ∇²ⁿ_2D::Vector{NF} = zeros(NF, trunc+2)              # initialized with zeros, ones
    ∇²ⁿ_2D_implicit::Vector{NF} = ones(NF, trunc+2)      # as this corresponds to no diffusion

    # PrimitiveEquation models need something more adaptive
    # and for each layer (to allow for varying orders/strength in the vertical)
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
    # always initialize the 2D arrays
    initialize!(scheme, model.time_stepping)
    
    # and the 3D arrays (different diffusion per layer) for primitive eq
    for k in 1:scheme.nlev 
        initialize!(scheme, k, model.geometry, model.time_stepping)
    end
end

"""$(TYPEDSIGNATURES)
Precomputes the 2D hyper diffusion terms in `scheme` based on the
model time step."""
function initialize!(   scheme::HyperDiffusion,
                        L::AbstractTimeStepper)

    (; trunc, ∇²ⁿ_2D, ∇²ⁿ_2D_implicit, power) = scheme
    (; Δt, radius) = L

    # time scale times 1/radius because time step Δt is scaled with 1/radius
    time_scale = scheme.time_scale.value/radius

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -trunc*(trunc+1)
    
    @inbounds for l in 0:trunc+1   # diffusion for every degree l, but indendent of order m
        eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

        # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
        ∇²ⁿ_2D[l+1] = -eigenvalue_norm^power/time_scale
        
        # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
        ∇²ⁿ_2D_implicit[l+1] = 1/(1-2Δt*∇²ⁿ_2D[l+1])           
    end
end

"""$(TYPEDSIGNATURES)
Precomputes the hyper diffusion terms in `scheme` for layer `k` based on the
model time step in `L`, the vertical level sigma level in `G`, and
the current (absolute) vorticity maximum level `vor_max`"""
function initialize!(   
    scheme::HyperDiffusion,
    k::Int,
    G::AbstractGeometry,
    L::AbstractTimeStepper,
    vor_max::Real = 0,
)
    (; trunc, resolution_scaling, ∇²ⁿ, ∇²ⁿ_implicit) = scheme
    (; power, power_stratosphere, tapering_σ) = scheme
    (; Δt, radius) = L
    σ = G.σ_levels_full[k]

    # Reduce diffusion time scale (=increase diffusion) with resolution
    # times 1/radius because time step Δt is scaled with 1/radius
    # time scale*3600 for [hrs] → [s]
    time_scale = scheme.time_scale.value/radius * (32/(trunc+1))^resolution_scaling

    # ADAPTIVE/FLOW AWARE
    # increase diffusion if maximum vorticity per layer is larger than scheme.vor_max
    if scheme.adaptive
        # /= as 1/time_scale*∇²ⁿ below
        time_scale /= 1 + (scheme.adaptive_strength-1)*max(0, vor_max/scheme.vor_max - 1)
    end

    # NORMALISATION
    # Diffusion is applied by multiplication of the eigenvalues of the Laplacian -l*(l+1)
    # normalise by the largest eigenvalue -lmax*(lmax+1) such that the highest wavenumber lmax
    # is dampened to 0 at the given time scale raise to a power of the Laplacian for hyperdiffusion
    # (=more scale-selective for smaller wavenumbers)
    largest_eigenvalue = -trunc*(trunc+1)
    
    # VERTICAL TAPERING for the stratosphere
    # go from 1 to 0 between σ=0 and tapering_σ
    tapering = max(0, (tapering_σ-σ)/tapering_σ)         # ∈ [0, 1]
    p = power + tapering*(power_stratosphere - power) 
        
    @inbounds for l in 0:trunc+1   # diffusion for every degree l, but indendent of order m
        eigenvalue_norm = -l*(l+1)/largest_eigenvalue   # normalised diffusion ∇², power=1

        # Explicit part (=-ν∇²ⁿ), time scales to damping frequencies [1/s] times norm. eigenvalue
        ∇²ⁿ[k][l+1] = -eigenvalue_norm^p/time_scale
        
        # and implicit part of the diffusion (= 1/(1-2Δtν∇²ⁿ))
        ∇²ⁿ_implicit[k][l+1] = 1/(1-2Δt*∇²ⁿ[k][l+1])           
    end
end

"""$(TYPEDSIGNATURES)
Pre-function to other `initialize!(::HyperDiffusion)` initialisors that
calculates the (absolute) vorticity maximum for the layer of `diagn`."""
function initialize!(   
    scheme::HyperDiffusion,
    diagn::DiagnosticVariablesLayer,
    G::AbstractGeometry,
    L::AbstractTimeStepper,
)
    scheme.adaptive || return nothing
    vor_min, vor_max = extrema(diagn.grid_variables.vor_grid)
    vor_abs_max = max(abs(vor_min), abs(vor_max))/G.radius
    initialize!(scheme, diagn.k, G, L, vor_abs_max)
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
                                ) where {NF<:AbstractFloat}
    lmax, mmax = matrix_size(tendency)      # 1-based
    @boundscheck size(tendency) == size(A) || throw(BoundsError)
    @boundscheck lmax <= length(∇²ⁿ_expl) == length(∇²ⁿ_impl) || throw(BoundsError)

    lm = 0
    @inbounds for m in 1:mmax   # loops over all columns/order m
        for l in m:lmax-1       # but skips the lmax+2 degree (1-based)
            lm += 1             # single index lm corresponding to harmonic l, m
            tendency[lm] = (tendency[lm] + ∇²ⁿ_expl[l]*A[lm])*∇²ⁿ_impl[l]
        end
        lm += 1             # skip last row for scalar quantities
    end
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity in the Barotropic models."""
function horizontal_diffusion!( diagn::DiagnosticVariablesLayer,
                                progn::PrognosticLayerTimesteps,
                                model::Barotropic,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    
    HD = model.horizontal_diffusion
    ∇²ⁿ = HD.∇²ⁿ_2D
    ∇²ⁿ_implicit = HD.∇²ⁿ_2D_implicit

    # Barotropic model diffuses vorticity (only variable)
    (; vor) = progn.timesteps[lf]
    (; vor_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion to vorticity and diffusion in the ShallowWater models."""
function horizontal_diffusion!( progn::PrognosticLayerTimesteps,
                                diagn::DiagnosticVariablesLayer,
                                model::ShallowWater,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    
    HD = model.horizontal_diffusion
    ∇²ⁿ = HD.∇²ⁿ_2D
    ∇²ⁿ_implicit = HD.∇²ⁿ_2D_implicit

    # ShallowWater model diffuses vorticity and divergence
    (; vor, div) = progn.timesteps[lf]
    (; vor_tend, div_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(div_tend, div, ∇²ⁿ, ∇²ⁿ_implicit)
end

"""$(TYPEDSIGNATURES)
Apply horizontal diffusion applied to vorticity, diffusion and temperature
in the PrimitiveEquation models. Uses the constant diffusion for temperature
but possibly adaptive diffusion for vorticity and divergence."""
function horizontal_diffusion!( progn::PrognosticLayerTimesteps,
                                diagn::DiagnosticVariablesLayer,
                                model::PrimitiveEquation,
                                lf::Int=1)      # leapfrog index used (2 is unstable)
    
    HD = model.horizontal_diffusion
    initialize!(HD, diagn, model.geometry, model.time_stepping)
    k = diagn.k                                 # current layer k
    ∇²ⁿ = HD.∇²ⁿ[k]                             # now pick operators at k
    ∇²ⁿ_implicit = HD.∇²ⁿ_implicit[k]

    # Primitive equation models diffuse vor and divergence more selective/adaptive
    (; vor, div, temp, humid) = progn.timesteps[lf]
    (; vor_tend, div_tend, temp_tend, humid_tend) = diagn.tendencies
    horizontal_diffusion!(vor_tend, vor, ∇²ⁿ, ∇²ⁿ_implicit)
    horizontal_diffusion!(div_tend, div, ∇²ⁿ, ∇²ⁿ_implicit)

    # but use the weaker normal diffusion for temperature, humidity
    ∇²ⁿ = HD.∇²ⁿ_2D
    ∇²ⁿ_implicit = HD.∇²ⁿ_2D_implicit
    horizontal_diffusion!(temp_tend, temp, ∇²ⁿ, ∇²ⁿ_implicit)
    model isa PrimitiveWet && horizontal_diffusion!(humid_tend, humid, ∇²ⁿ, ∇²ⁿ_implicit)
end